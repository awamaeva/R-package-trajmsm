#' @title Inverse Probability Weighting
#' @description Compute stabilized and unstabilized weights, with or without censoring.
#' @name inverse_probability_weighting
#' @param numerator To choose between stabilized and unstabilized weights.
#' @param identifier Name of the column of the unique identifier.
#' @param baseline Name of the baseline covariates.
#' @param covariates Name of the time-varying covariates.
#' @param treatment Name of the time-varying treatment.
#' @param include_censor Logical value TRUE/FALSE to include or not a censoring variable.
#' @param censor Name of the censoring variable.
#' @param obsdata Observed data in wide format.
#' @return Inverse Probability Weights (Stabilized and Unstabilized) with and without censoring.
#' @export
#' @author Awa Diop, Denis Talbot
#' @examples
#' obsdata = gendata(n = 1000, format = "wide",total_followup = 3, seed = 945)
#' baseline_var <- c("age","sex")
#' covariates <- list(c("hyper2011", "bmi2011"),
#' c("hyper2012", "bmi2012"),c("hyper2013", "bmi2013"))
#' treatment_var <- c("statins2011","statins2012","statins2013")
#' stabilized_weights = inverse_probability_weighting(numerator = "stabilized",
#' identifier = "id", covariates = covariates, treatment = treatment_var,
#' baseline = baseline_var, obsdata = obsdata)

inverse_probability_weighting <- function(numerator = c("stabilized", "unstabilized"), identifier,
                                          baseline, covariates, treatment,
                                           include_censor = FALSE, censor,
                                          obsdata){
  if (!is.data.frame(obsdata)) {
    stop("obsdata must be a data frame in wide format")
  }

  # Validate numerator
  numerator <- match.arg(numerator)

  ipw_result <- switch(numerator,
                       "unstabilized" = ifelse(include_censor,
                                               unstabilized_ipcw(identifier = identifier, treatment = treatment,
                                                                 covariates = covariates, baseline = baseline,
                                                                 censor = censor, obsdata = obsdata),
                                               unstabilized_iptw(identifier = identifier, treatment = treatment,
                                                                 covariates = covariates, baseline = baseline,  obsdata = obsdata)),
                       "stabilized" = ifelse(include_censor,
                                             stabilized_ipcw(identifier = identifier, treatment = treatment,
                                                             covariates = covariates, baseline = baseline, censor = censor,  obsdata = obsdata),
                                             stabilized_iptw(identifier = identifier, treatment = treatment, covariates = covariates,
                                                             baseline = baseline,  obsdata = obsdata))
  )

  return(ipw_result)
}
