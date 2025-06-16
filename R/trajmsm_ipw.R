#' @title Marginal Structural Model and Latent Class of Growth Analysis estimated with IPW
#' @description Estimate parameters of LCGA-MSM using IPW.
#' @name trajmsm_ipw
#' @param formula1 Specification of the model for the outcome to be fitted for a binomial or gaussian distribution.
#' @param formula2 Specification of the model for the outcome to be fitted for a survival outcome.
#' @param family Specification of the error distribution and link function to be used in the model.
#' @param identifier Name of the column of the unique identifier.
#' @param treatment Time-varying treatment.
#' @param baseline Name of the baseline covariates.
#' @param covariates Names of the time-varying covariates (should be a list).
#' @param obsdata Dataset to be used in the analysis.
#' @param numerator Type of weighting ("stabilized" or "unstabilized").
#' @param weights A vector of estimated weights. If NULL, the weights are computed by the function \code{IPW}.
#' @param censor Name of the censoring variable.
#' @param include_censor Logical, if TRUE, includes censoring.
#' @param treshold For weight truncation.
#' @return Provides a matrix of estimates for LCGA-MSM, obtained using IPW.
#' @export trajmsm_ipw
#' @importFrom survival coxph
#' @import flexmix sandwich
#' @importFrom stats na.omit rbinom plogis qlogis  reshape glm
#' binomial coef as.formula ave aggregate relevel pnorm sd quantile model.matrix
#' @return Provides a matrix of estimates for LCGA-MSM, obtained using IPW.
#' @examples
#' \donttest{
#' obsdata_long = gendata(n = 1000, format = "long", total_followup = 6, seed = 845)
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
#'            treatment = treatment_var, family = "binomial",
#'            obsdata = obsdata,numerator = "stabilized", include_censor = FALSE, treshold = 0.99)
#'resmsm_ipw
#'}

trajmsm_ipw <- function(formula1, formula2, family, identifier, treatment, covariates,
                        baseline, obsdata,
                        numerator = "stabilized", include_censor = FALSE,
                        censor, weights = NULL, treshold = 0.99) {

  # Compute weights if not provided
  if (is.null(weights)) {
    weights <- inverse_probability_weighting(identifier = as.character(identifier), covariates = covariates,
                                             treatment = treatment, baseline = baseline,
                                             numerator = numerator,

                                             include_censor = include_censor, censor = censor, obsdata = obsdata)[[1]][, length(treatment)]

    obsdata$weights <- ifelse(quantile(weights, treshold, na.rm = TRUE) > weights, quantile(weights, treshold, na.rm = TRUE), weights)
  } else {
    obsdata$weights <- weights
  }

  # Model fitting
  if (family == "survival") {
    # Prepare the list of arguments for coxph
    args_list <- list(
      formula = formula2,
      data = obsdata,
      cluster = as.formula(paste0("~", identifier)), # Dynamically construct cluster formula
      weights = weights
    )

    # Use do.call to invoke coxph with dynamically constructed arguments
    mod_glm <- do.call("coxph", args_list)
  } else {
    # Fit a standard glm model
    mod_glm <- glm(formula = formula1, data = obsdata, weights = obsdata$weights, family = family)

    # Compute robust standard errors using the sandwich package
    cluster_formula <- as.formula(paste0("~", identifier))
    robust_se <- sqrt(diag(vcovCL(mod_glm, cluster = cluster_formula)));
  }

  # Extracting model results
  coefs <- coef(mod_glm)
  if (family == "survival") {
    # Extract standard errors for survival models
    se <- coef(summary(mod_glm))[, 2]
  } else {
    se <- robust_se  # Use robust standard errors for non-survival models
  }
  pvalue <- 2 * pnorm(-abs(coefs) / se)

  ic_lo <- coefs - 1.96 * se
  ic_up <- coefs + 1.96 * se

  restrajmsm_ipw <- cbind(coefs, se, pvalue, ic_lo, ic_up)
  colnames(restrajmsm_ipw) <- c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI")

  return(restrajmsm_ipw)
}
