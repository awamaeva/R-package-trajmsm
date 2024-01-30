#' @title Stabilized Inverse Probability of Censoring Weights
#' @description Compute stabilized Inverse Probability of Censoring Weights
#' @name stabilized_ipcw
#' @param identifier Name of the column for unique identifier.
#' @param treatment Time-varying treatment.
#' @param covariates List of time-varying covariates.
#' @param baseline Baseline covariates.
#' @param total_follow_up Total length of follow-up.
#' @param censor Name of the censoring variable.
#' @param obsdata Observed data in wide format.
#' @return Stabilized Inverse Probability of Censoring Weights
#' @noRd
#' @export
#' @author Awa Diop, Denis Talbot
#' @note This function requires data in a wide format.

stabilized_ipcw <- function(identifier, treatment, covariates, baseline, censor, total_follow_up, obsdata) {
  if (!is.data.frame(obsdata)) {
    stop("obsdata must be a data frame in wide format")
  }

  required_args <- c("identifier", "baseline", "covariates", "treatment", "total_follow_up", "censor", "obsdata")
  missing_args <- setdiff(required_args, names(match.call()))
  if (length(missing_args) > 0) {
    stop(paste(missing_args, collapse = ", "), " not specified")
  }

  sw_treatment_temp <- matrix(numeric(0), nrow = nrow(obsdata), ncol = length(treatment))
  sw_censor_temp <- matrix(numeric(0), nrow = nrow(obsdata), ncol = length(treatment))

  for (i in 2:length(treatment)) {
    past_treatments <- paste0(treatment[1:(i-1)], collapse = "+")
    current_covariates <- paste0(unlist(covariates[1:i]), collapse = "+")
    baseline_covariates <- paste0(baseline, collapse = "+")

    form_num_treatment <- as.formula(paste(paste0(treatment[i], "~", past_treatments)))
    form_num_censor <- as.formula(paste(paste0(censor[i], "~ 1")))

    form_denom_treatment <- as.formula(paste(paste0(treatment[i], "~", past_treatments, "+", current_covariates, "+", baseline_covariates)))
    form_denom_censor <- as.formula(paste(paste0(censor[i], "~", past_treatments, "+", current_covariates, "+", baseline_covariates)))

    fit_num_treatment <- glm(form_num_treatment, family = binomial(link = "logit"), data = obsdata)
    fit_num_censor <- glm(form_num_censor, family = binomial(link = "logit"), data = obsdata)

    fit_denom_treatment <- glm(form_denom_treatment, family = binomial(link = "logit"), data = obsdata)
    fit_denom_censor <- glm(form_denom_censor, family = binomial(link = "logit"), data = obsdata)

    ps_treatment <- predict(fit_num_treatment, type = "response")
    ps_censor <- predict(fit_num_censor, type = "response")


    sw_treatment_temp[, i][obsdata[,censor[i-1]] == 0] <-  ((1 - obsdata[, treatment[i]][obsdata[,censor[i-1]] == 0]) * (1 - ps_treatment[obsdata[,censor[i-1]] == 0])) / (1 - predict(fit_denom_treatment, type = "response")) +
      (obsdata[, treatment[i]][obsdata[,censor[i-1]] == 0] * ps_treatment[obsdata[,censor[i-1]] == 0]) / predict(fit_denom_treatment, type = "response")

    sw_censor_temp[, i][obsdata[,censor[i-1]] == 0] <- ((1 - ps_censor[obsdata[,censor[i-1]] == 0]) / predict(fit_denom_censor, type = "response"))
  }

  # Calculate stabilized weights for the first time point
  # Formulating numerator model for treatment and censoring at t=1
  form_num_treatment_t1 <- as.formula(paste0(treatment[1], "~ 1"))
  form_num_censor_t1 <- as.formula(paste0(censor[1], "~ 1"))

  # Formulating denominator model for treatment and censoring at t=1
  # Including only the baseline covariates since it's the first time point
  form_denom_treatment_t1 <- as.formula(paste0(treatment[1], "~", paste0(baseline, collapse = "+")))
  form_denom_censor_t1 <- as.formula(paste0(censor[1], "~", paste0(baseline, collapse = "+")))

  # Fitting models for numerator at t=1
  fit_num_treatment_t1 <- glm(form_num_treatment_t1, family = binomial(link = "logit"), data = obsdata)
  fit_num_censor_t1 <- glm(form_num_censor_t1, family = binomial(link = "logit"), data = obsdata)

  # Fitting models for denominator at t=1
  fit_denom_treatment_t1 <- glm(form_denom_treatment_t1, family = binomial(link = "logit"), data = obsdata)
  fit_denom_censor_t1 <- glm(form_denom_censor_t1, family = binomial(link = "logit"), data = obsdata)

  # Predicting probabilities for treatment and censoring at t=1
  ps_treatment_t1 <- predict(fit_num_treatment_t1, type = "response")
  ps_censor_t1 <- predict(fit_num_censor_t1, type = "response")

  # Calculating stabilized weights for the first time point
  sw_treatment_temp[, 1] <- (obsdata[, treatment[1]] * ps_treatment_t1) / predict(fit_denom_treatment_t1, type = "response")
  sw_censor_temp[, 1] <- (1 - ps_censor_t1) / predict(fit_denom_censor_t1, type = "response")

  # Compute final stabilized weights
  sw_censor_treatment <- sw_treatment_temp * sw_censor_temp

  weights <- t(apply(sw_treatment_temp * sw_censor_temp, 1, cumprod))

  return(list(weights))
}
