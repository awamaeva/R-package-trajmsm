#' @title Marginal Structural Model and Latent Class of Growth Analysis estimated with IPW
#' @description Combine Marginal Structural Model and Latent Class of Growth Analysis
#' @name trajmsm_ipw
#' @param formula1 Specification of the model for the outcome to be fitted for a binomial or gaussian distribution.
#' @param formula2 Specification of the model for the outcome to be fitted for a survival outcome.
#' @param family Specification of the error distribution and link function to be used in the model.
#' @param identifier Name of the column for unique identification.
#' @param treatment Time-varying treatment.
#' @param baseline Baseline covariates.
#' @param covariates Time-varying covariates.
#' @param total_followup Number of measuring times.
#' @param number_traj An integer to fix the number of trajectory groups.
#' @param obsdata Dataset to be used in the analysis.
#' @param numerator Type of weighting ("stabilized" or "unstabilized").
#' @param weights A vector of estimated weights. If NULL, the weights are computed by the function \code{IPW}.
#' @param treshold For weight truncation.
#' @return Stabilized and unstabilized inverse of probabilities
#' @export trajmsm_ipw
#' @import sandwich
#' @importFrom survival coxph
#' @import flexmix
#' @examples
#' obsdata_long = gendata(n = 1000, format = "long", total_followup = 6, seed = 945)
#' years <- 2011:2016
#' baseline_var <- c("age","sex")
#' variables <- c("hyper", "bmi")
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
#' trajmsm_long$ipw_group <- relevel(trajmsm_long$class, ref = "1")
#' obsdata = reshape(data = trajmsm_long, direction = "wide", idvar = "id",
#' v.names = c("statins","bmi","hyper"), timevar = "time", sep ="")
#' formula = paste0("y ~", paste0(treatment_var,collapse = "+"), "+",
#'                 paste0(unlist(covariates), collapse = "+"),"+",
#'                 paste0(baseline_var, collapse = "+"))
#'
#'resmsm_ipw = trajmsm_ipw(formula1 = as.formula("y ~ ipw_group"),
#'            identifier = "id", baseline = baseline_var, covariates = covariates,
#'            treatment = treatment_var, number_traj=3,total_followup = 6, family = "binomial",
#'            obsdata = obsdata,numerator = "stabilized", include_censor = FALSE)
#'resmsm_ipw


trajmsm_ipw <- function(formula1, formula2, family, identifier, treatment, covariates,
                         baseline, total_followup, number_traj, obsdata,
                         numerator = "stabilized",include_censor, censor, weights = NULL, treshold = 0.99) {

  # Compute weights if not provided
  if (is.null(weights)) {
    weights <- inverse_probability_weighting(identifier = identifier, covariates = covariates,
                                             treatment = treatment, baseline = baseline,
                                             total_followup = total_followup, numerator = numerator,
                                             include_censor = include_censor, censor = censor,obsdata = obsdata)[[1]][, total_followup]

    obsdata$weights <- ifelse(quantile(weights, treshold, na.rm = TRUE)> weights, quantile(weights, treshold, na.rm = TRUE), weights)
  } else {
    obsdata$weights <- weights
  }

  # Model fitting
  cluster_formula <- as.formula(paste0("~", identifier))
  if (family == "gaussian") {
    mod_glm <- glm(formula1, data = obsdata, weights = weights, family = gaussian)
  } else if (family == "binomial") {
    mod_glm <- glm(formula1, data = obsdata, weights = weights, family = binomial)
  } else if (family == "survival") {
    mod_glm <- coxph(formula2, data = obsdata,cluster = get(identifier), weights = weights)
  }

  # Extracting model results
  coefs <- coef(summary(mod_glm))[, 1]
  se <- sqrt(diag(vcovCL(mod_glm, cluster = cluster_formula)))
  pvalue <- coef(summary(mod_glm))[,4]
  ic_lo <- coefs - 1.96 * se
  ic_up <- coefs + 1.96 * se

  res_trajmsm <- cbind(coefs, se, pvalue, ic_lo, ic_up)
  colnames(res_trajmsm) <- c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI")

  return(res_trajmsm)
}


