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
#' @param y Outcome of interest.
#' @param total_followup Number of measuring times.
#' @param number_traj An integer to fix the number of trajectory groups.
#' @param obsdata Dataset to be used in the analysis.
#' @param numerator Type of weighting ("stabilized" or "unstabilized").
#' @param weights A vector of estimated weights. If NULL, the weights are computed by the function \code{IPW}.
#' @param treshold For weight truncation.
#' @return Stabilized and unstabilized inverse of probabilities
#' @export trajmsm_ipw
#' @importFrom sandwich
#' @importFrom survival coxph
#' @import flexmix
#' @examples
# Example usage of the function
#' Obsdata_long = gendata_trajmsm(n = 1000, include_censor = TRUE, format = "long", seed = 945)
#' baseline_var <- c("age","sex")
#' covariates <- list(c("hyper2011", "bmi2011"),c("hyper2012", "bmi2012"),c("hyper2013", "bmi2013"))
#' treatment_var <- c("statins2011","statins2012","statins2013")
#' censor_var = c("censor2011","censor2012","censor2013")
#' formula = as.formula(cbind(statins, 1 - statins) ~ time)
#' restraj = build_traj(obsdata = Obsdata_long, number_traj = 3, formula = formula, identifier = "id")
#' Datapost = restraj$data_post
#' trajmsm_long <- merge(Obsdata_long, Datapost, by = "id")
#'     AggFormula <- as.formula(paste("statins", "~", "time", "+", "class"))
#'     AggTrajData <- aggregate(AggFormula, data = trajmsm_long, FUN = mean)
#'     AggTrajData
#'trajmsm_long[ , "traj_group"] <- factor(ifelse(trajmsm_long[ , "class"] == "2" ,"Group1" ,
#' ifelse (trajmsm_long[ , "class"]== "1" , "Group2" ,"Group3")))
#' trajmsm_long[ , "traj_group"] <- relevel(trajmsm_long[ , "traj_group"], ref = "Group3")
#' trajmsm_wide = reshape(trajmsm_long, direction = "wide", idvar = "id",
#' v.names = c("statins","bmi","hyper","censor"), timevar = "time", sep ="")
#'trajmsm_ipw(formula1 = as.formula("y ~ traj_group"),
#'            identifier = "id", baseline = baseline_var, covariates = covar, treatment = treatment_var,
#'            y="y", number_traj=3,total_followup = 3, family = "binomial",
#'            obsdata = trajmsm_wide,numerator = "stabilized", censor = censor_var, include_censor = TRUE)


trajmsm_ipw <- function(formula1, formula2, family, identifier, treatment, covariates,
                         baseline, y, total_followup, number_traj, obsdata,
                         numerator = "stabilized",include_censor, censor, weights = NULL, treshold = 0.99) {

  # Compute weights if not provided
  if (is.null(weights)) {
    weights <- inverse_probability_weighting(identifier = identifier, covariates = covariates,
                                             treatment = treatment, baseline = baseline,
                                             total_follow_up = total_followup, numerator = numerator,
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


