#' @title Marginal Structural Model and Latent Class of Growth Analysis estimated with IPW
#' @description Combine Marginal Structural Model and Latent Class of Growth Analysis
#' @name trajhrmsm_ipw
#' @param formula specification of the model for the outcome to be fitted for a binomial or gaussian distribution.
#' @param family specification of the error distribution and link function to be used in the model.
#' @param degree_traj to specify the polynomial degree for modelling the time-varying treatment.
#' @param identifier  name of the column for unique identifiant.
#' @param baseline name of baseline covariates.
#' @param covariates names of time-varying covariates in a wide format.
#' @param treatment name of time-varying treatment.
#' @param var_cov names of the time-varying covariates in a long format.
#' @param ntimes_interval length of a time-interval (s).
#' @param total_followup total length of follow-up.
#' @param censor name of the censoring variable.
#' @param number_traj number of trajectory groups.
#' @param weights a vector of estimated weights. If NULL, the weights are computed by the function.
#' @param obsdata data in a long format.
#' @author Awa Diop, Denis Talbot
#' @export trajhrmsm_ipw
#' @importFrom sandwich
#' @import flexmix
#' @examples
#' Obsdata_long = gendata(n = 1000, format = "long", total_followup = 5, seed = 945)
#' formula = as.formula("y ~ factor(traj) + factor(Interval)")
#' baseline_var <- c("age","sex")
#' covariates <- list(c("hyper2011", "bmi2011"),c("hyper2012", "bmi2012"),c("hyper2013", "bmi2013"))
#' treatment_var <- c("statins2011","statins2012","statins2013")
#' var_cov <- c("statins","hyper", "bmi")
#' resipw <- trajhrmsm_ipw(formula = formula, degree_traj = "linear", numerator = "stabilized",
#' identifier = "id", baseline = baseline_var, covariates = covariates, treatment = treatment_var,
#' outcome = "y", name_traj = "traj", name_interv = "Interval",var_cov= c("statins","hyper","bmi"),
#' include_censor = FALSE, ntimes_interval = 3,total_followup = 5, tim = "time",family = "poisson",
#' number_traj = 3, obsdata = Obsdata_long, treshold = 0.999)
#' resipw$res_trajHRMSM_IPW



trajhrmsm_ipw <- function(formula ,
                          degree_traj = c("linear","quadratic","cubic"),
                          numerator = c("stabilized", "unstabilized"),
                          identifier, baseline, covariates, treatment, outcome, name_traj, name_interv,var_cov,
                          include_censor = FALSE, ntimes_interval,total_followup,time,family = "poisson",censor = censor,
                          number_traj, obsdata, weights = NULL, treshold = 0.999){

  if(is.null(weights)){
    obsdata_wide <- reshape(obsdata, direction = "wide", idvar = identifier, v.names = var_cov, timevar = time, sep ="")
    # Compute the weights for all different regimes of treatment
    Weights = inverse_probability_weighting(identifier = identifier, covariates = covariates,
                                            treatment = treatment, baseline = baseline,
                                            total_follow_up = total_followup, numerator = numerator,
                                            include_censor = include_censor, censor = censor,obsdata = obsdata_wide)[[1]];

      weights_trunc <- sapply(1:ntimes_interval, function(x){
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
  #Choice of degree for the polynomial form to buil the trajectory groups
  if(degree_traj == "linear"){

  res_traj = build_traj(obsdata  = na.omit(dat_sub[,c(treatment_name,"time2","identifier2")]),
                        number_traj = number_traj,formula = as.formula(paste("cbind(", treatment_name, ", 1 -", treatment_name, ") ~ time2")),
                        identifier = "identifier2")
  }

  if(degree_traj == "quadratic"){

    res_traj = build_traj(obsdata = na.omit(dat_sub[,c(treatment_name,"time2","identifier2")]), number_traj = number_traj,
                          formula  = as.formula(paste("cbind(", treatment_name, ", 1 -", treatment_name, ") ~ time2 + I(time2^2)")),
                          identifier = "identifier2")
  }

  if(degree_traj == "cubic"){

    res_traj = build_traj(obsdata = na.omit(dat_sub[,c(treatment_name,"time2","identifier2")]), number_traj = number_traj,
                          formula  = as.formula(paste("cbind(", treatment_name, ", 1 -", treatment_name, ") ~ time2 + I(time2^2) + I(time2^3)")),
                          identifier = "identifier2")
  }

  if(length(unique(res_traj$data_post$class)) < number_traj){stop("number of trajectory groups identified is inferior to the target number.")}

  if(length(unique(res_traj$data_post$class)) == number_traj){
    dclass <- data.frame(traj = factor(res_traj$data_post[,"class"]), identifier2 = res_traj$data_post[,"identifier2"])
    colnames(dclass)[1] <- name_traj
    dat_final <- merge(dat_sub, dclass, by = "identifier2")

    cluster_formula <- as.formula(paste0("~", "identifier2"))

    mean_adh <- aggregate(as.formula(paste0(treatment_name, "~", name_traj)), FUN = mean, data = dat_final)
    ord_adh<- order(-mean_adh[,2])
    ref <- ord_adh[length(ord_adh)]
    dat_final[ ,name_traj] <- relevel(factor(dat_final[ ,name_traj]), ref = ref)
    colnames(dat_final)[which(names(dat_final) == "Interv")]<- name_interv

    mod_glm = glm(formula = formula, weights = IPW,family = family,data=dat_final);
    coefs <- summary(mod_glm)$coefficients[1:number_traj,1];
    se   <- sqrt(diag(vcovCL(mod_glm, cluster = cluster_formula)))[1:number_traj];
    pvalue <- round(coef(summary(mod_glm))[1:number_traj,4],4)
    IClo = coefs- 1.96*se ;
    ICup = coefs + 1.96*se;

    res_trajHRMSM_IPW = cbind(coefs,se,   pvalue, IClo, ICup);
    colnames(res_trajHRMSM_IPW) = c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI");
    res_trajHRMSM_IPW = res_trajHRMSM_IPW[2:number_traj,]

return(list(res_trajHRMSM_IPW = res_trajHRMSM_IPW, res_traj = res_traj, mean_adh = mean_adh))}
}
