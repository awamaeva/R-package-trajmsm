#' @title Parametric g-formula
#' @description Estimate parameters of LCGA-MSM using g-formula
#'  and bootstrap to get standard errors.
#' @name trajmsm_gform
#' @param formula Specification of the model for the outcome to be fitted.
#' @param identifier Name of the column of the unique identifier.
#' @param baseline  Vector of names of the baseline covariates.
#' @param covariates List of names of the time-varying covariates.
#' @param treatment Vector of names of the time-varying treatment.
#' @param outcome Name of the outcome of interest.
#' @param total_followup Total length of follow-up.
#' @param time Name of the time variable.
#' @param time_values Measuring times.
#' @param rep Number of repetitions for the bootstrap.
#' @param trajmodel Trajectory model built with the observed treatment.
#' @param var_cov Names of the time-varying covariates.
#' @param ref The reference trajectory group.
#' @param obsdata Observed data in wide format.
#' @return Provides a matrix of estimates for LCGA-MSM, obtained using the g-formula method.
#' @importFrom stats quasibinomial var
#' @importFrom utils combn
#' @export
#' @examples
#' \donttest{
#' obsdata_long = gendata(n = 1000, format = "long", total_followup = 6, seed = 845)
#' years <- 2011:2016
#' baseline_var <- c("age","sex")
#' variables <- c("hyper", "bmi")
#' var_cov <- c("statins","hyper", "bmi")
#' covariates <- lapply(years, function(year) {
#' paste0(variables, year)})
#' treatment_var <- paste0("statins", 2011:2016)
#' formula_treatment = as.formula(cbind(statins, 1 - statins) ~ time)
#' restraj = build_traj(obsdata = obsdata_long, number_traj = 3,
#' formula = formula_treatment, identifier = "id")
#' datapost = restraj$data_post
#' trajmsm_long <- merge(obsdata_long, datapost, by = "id")
#'     AggFormula <- as.formula(paste("statins", "~", "time", "+", "class"))
#'     AggTrajData <- aggregate(AggFormula, data = trajmsm_long, FUN = mean)
#'     AggTrajData
#' obsdata = reshape(data = trajmsm_long, direction = "wide", idvar = "id",
#' v.names = c("statins","bmi","hyper"), timevar = "time", sep ="")
#'formula = paste0("y ~", paste0(treatment_var,collapse = "+"), "+",
#'                 paste0(unlist(covariates), collapse = "+"),"+",
#'                 paste0(baseline_var, collapse = "+"))
#'resmsm_gform <- trajmsm_gform(formula = formula, identifier = "id",rep = 5,
#'baseline = baseline_var, covariates = covariates, var_cov = var_cov,
#'treatment = treatment_var, outcome = "y", total_followup = 6,time = "time",
#' time_values = years, trajmodel = restraj$traj_model,ref = "1", obsdata =   obsdata )
#' resmsm_gform
#' }
#' @author Awa Diop Denis Talbot

trajmsm_gform <- function(formula = formula, rep = 50,
                                  identifier,baseline,covariates,treatment,outcome,
                                  total_followup,time = time,time_values, var_cov,
                                  trajmodel,ref, obsdata){


  stopifnot(!is.null(identifier));
  stopifnot(!is.null(baseline));
  stopifnot(!is.null(covariates));
  stopifnot(!is.null(treatment));
  stopifnot(!is.null(outcome));
  stopifnot(!is.null(total_followup));
  stopifnot(!is.null(rep));
  stopifnot(!is.null(obsdata));
  stopifnot(!is.null(trajmodel));
  stopifnot(!is.null(time));

    df <- data.frame()
    bootf=function(df,x=index){
      #Echantillons bootstrap
      df=obsdata[x,];

      res = gformula(formula = formula,outcome = outcome, treatment = treatment, covariates = covariates,
                     baseline = baseline,ntimes_interval = total_followup, obsdata = df)$counter_means
      colnames(res) <- "Y";

      #Counterfactual means + trajectory groups
      obsdataG = data.frame(res);
      treatment_names <- sub("\\d+", "", treatment)
      treatment_name <- unique(treatment_names)[1]
      obsdataG$gform_group = factor(predict_traj(identifier = identifier, total_followup = total_followup,
                                                 treatment = treatment_name, time = time,time_values = time_values,
                                                 trajmodel = trajmodel)$post_class);

      obsdataG$gform_group <- relevel(as.factor(obsdataG$gform_group), ref = ref)
      #Estimation
      mod = summary(glm(Y ~ gform_group, family = quasibinomial, data = obsdataG));
      return(coef(mod)[,1]);
    }

    #Bootstrap
    list_res <- list()
    for(b in 1:rep){
      index <- sample(1:nrow(obsdata), nrow(obsdata) ,replace = T)
      list_res[b] <- list(bootf(df,x = index))
    }

    result.coef.boot <- do.call(rbind,list_res)
    #results
    res = gformula(formula = formula,outcome = outcome, treatment = treatment, covariates = covariates,baseline = baseline,
                    ntimes_interval = total_followup, obsdata = obsdata)$counter_means
    colnames(res) <- "Y";

    #Counterfactual means + trajectory groups
    obsdataG = data.frame(res);
    treatment_names <- sub("\\d+", "", treatment)
    treatment_name <- unique(treatment_names)[1]
    obsdataG$gform_group = factor(predict_traj(identifier = identifier, total_followup = total_followup,
                                               treatment = treatment_name, time = time,time_values = time_values,
                                               trajmodel = trajmodel)$post_class);
    obsdataG$gform_group <- relevel(as.factor(obsdataG$gform_group), ref = ref)
    #Estimation
    mod = summary(glm(Y ~ gform_group, family = quasibinomial, data = obsdataG));
    coefs = mod$coefficients[,1]
    se = apply(result.coef.boot,2,sd)
    pvalue <- 2*pnorm(-abs(coefs)/se)
    lo.ci = coefs - 1.96*se
    up.ci = coefs + 1.96*se
    results_msm_gform = rbind(coefs,se,pvalue,lo.ci, up.ci);
    rownames(results_msm_gform) = c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI")
    return(t(results_msm_gform));
}
