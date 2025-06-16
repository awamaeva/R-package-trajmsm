#' @title History Restricted MSM and Latent Class of Growth Analysis estimated with IPW.
#' @description Estimate parameters of LCGA-HRMSM using IPW.
#' @name trajhrmsm_ipw
#' @param family specification of the error distribution and link function to be used in the model.
#' @param numerator To choose between stabilized and unstabilized weights.
#' @param degree_traj To specify the polynomial degree for modelling the time-varying treatment.
#' @param identifier  Name of the column of the unique identifier.
#' @param baseline Names of the baseline covariates.
#' @param covariates Names of the time-varying covariates (should be a list).
#' @param treatment Name of the time-varying treatment.
#' @param outcome Name of the outcome variable.
#' @param var_cov Names of the time-varying variables.
#' @param ntimes_interval Length of a time-interval (s).
#' @param total_followup Total length of follow-up.
#' @param censor Name of the censoring variable.
#' @param include_censor Logical, if TRUE, includes censoring.
#' @param number_traj Number of trajectory groups.
#' @param weights A vector of estimated weights. If NULL, the weights are computed by the function.
#' @param time_values Values of the time variable.
#' @param time Name of the time variable.
#' @param treshold For weight truncation.
#' @param obsdata Data in a long format.
#' @return Provides a matrix of estimates for LCGA-HRMSM, obtained using IPW.
#' @author Awa Diop, Denis Talbot
#' @import sandwich
#' @import flexmix
#' @importFrom stats na.omit rbinom plogis qlogis  reshape glm
#' binomial coef as.formula ave aggregate relevel pnorm sd quantile model.matrix
#' @export trajhrmsm_ipw
#' @examples
#' \donttest{
#' obsdata_long = gendata(n = 5000, format = "long", total_followup = 8,
#' timedep_outcome = TRUE,  seed = 845)
#' baseline_var <- c("age","sex")
#' years <- 2011:2018
#' variables <- c("hyper", "bmi")
#' covariates <- lapply(years, function(year) {
#' paste0(variables, year)})
#' treatment_var <- paste0("statins", 2011:2018)
#' var_cov <- c("statins","hyper", "bmi","y")
#' reshrmsm_ipw <- trajhrmsm_ipw(degree_traj = "linear", numerator = "stabilized",
#' identifier = "id", baseline = baseline_var,
#' covariates = covariates, treatment = treatment_var,
#' outcome = "y", var_cov= var_cov,include_censor = FALSE,
#'  ntimes_interval = 6,total_followup = 8, time = "time", time_values = 2011:2018,
#' family = "poisson", number_traj = 3, obsdata = obsdata_long, treshold = 1)
#' reshrmsm_ipw$res_trajhrmsm_ipw
#' }



trajhrmsm_ipw <- function(degree_traj = c("linear","quadratic","cubic"),
                          numerator = c("stabilized", "unstabilized"),
                          identifier, baseline, covariates, treatment, outcome,var_cov,
                          include_censor = FALSE, ntimes_interval,total_followup,time,time_values, family = "poisson",censor = censor,
                          number_traj, obsdata, weights = NULL, treshold = 0.999){


  if(is.null(weights)){
    obsdata_wide <- reshape(obsdata, direction = "wide", idvar = identifier, v.names = var_cov, timevar = time, sep ="")
    # Compute the weights for all different regimes of treatment
    Weights = inverse_probability_weighting(identifier = identifier, covariates = covariates,
                                            treatment = treatment, baseline = baseline,
                                            numerator = numerator,
                                            include_censor = include_censor, censor = censor,obsdata = obsdata_wide)[[1]];

      weights_trunc <- sapply(1:total_followup, function(x){
      weights <- ifelse(quantile(Weights[, x], treshold, na.rm = TRUE)> Weights[, x], quantile(Weights[, x], treshold, na.rm = TRUE), Weights[, x])
      return(weights)
    })
    }

  if(!is.null(weights)){
    weights_trunc <- weights
  }


  obsdata$IPW <- as.vector(weights_trunc)

  dat_sub = data.frame(do.call(rbind, split_data(obsdata = obsdata, total_followup = total_followup,
                                                 ntimes_interval = ntimes_interval,
                                                 time = time, time_values = time_values, identifier = identifier)))

  treatment_names <- sub("\\d+", "", treatment)
  treatment_name <- unique(treatment_names)[1]
  #Choice of degree for the polynomial form to build the trajectory groups
  if(degree_traj == "linear"){

  restraj = build_traj(obsdata  = na.omit(dat_sub[,c(treatment_name,"time2","identifier2")]),
                        number_traj = number_traj,formula = as.formula(paste("cbind(", treatment_name, ", 1 -", treatment_name, ") ~ time2")),
                        identifier = "identifier2")
  }

  if(degree_traj == "quadratic"){

    restraj = build_traj(obsdata = na.omit(dat_sub[,c(treatment_name,"time2","identifier2")]), number_traj = number_traj,
                          formula  = as.formula(paste("cbind(", treatment_name, ", 1 -", treatment_name, ") ~ time2 + I(time2^2)")),
                          identifier = "identifier2")
  }

  if(degree_traj == "cubic"){

    restraj = build_traj(obsdata = na.omit(dat_sub[,c(treatment_name,"time2","identifier2")]), number_traj = number_traj,
                          formula  = as.formula(paste("cbind(", treatment_name, ", 1 -", treatment_name, ") ~ time2 + I(time2^2) + I(time2^3)")),
                          identifier = "identifier2")
  }

  if(length(unique(restraj$data_post$class)) < number_traj){stop("number of trajectory groups identified is inferior to the target number.")}

  if(length(unique(restraj$data_post$class)) == number_traj){
    dclass <- data.frame(ipw_group = factor(restraj$data_post[,"class"]), identifier2 = restraj$data_post[,"identifier2"])
    dat_final <- merge(dat_sub, dclass, by = "identifier2")

    cluster_formula <- as.formula(paste0("~", "identifier2"))

    mean_adh <- aggregate(as.formula(paste0(treatment_name, "~", "ipw_group")), FUN = mean, data = dat_final)
    ord_adh<- order(-mean_adh[,2])
    ref <- ord_adh[length(ord_adh)]
    dat_final[ ,"ipw_group"] <- relevel(factor(dat_final[ ,"ipw_group"]), ref = ref)


    mod_glm = glm(formula =as.formula(paste(outcome, "~ factor(ipw_group) + factor(Interv)")), weights = dat_final[,"IPW"],family = family,data=dat_final);
    coefs <- summary(mod_glm)$coefficients[1:number_traj,1];
    se   <- sqrt(diag(vcovCL(mod_glm, cluster = cluster_formula)))[1:number_traj];
    pvalue <- 2*pnorm(-abs(coefs)/se)
    IClo = coefs- 1.96*se ;
    ICup = coefs + 1.96*se;

    res_trajhrmsm_ipw = cbind(coefs,se,  pvalue, IClo, ICup);
    colnames(res_trajhrmsm_ipw) = c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI");


return(list( res_trajhrmsm_ipw =  res_trajhrmsm_ipw, restraj = restraj, mean_adh = mean_adh))}
}
