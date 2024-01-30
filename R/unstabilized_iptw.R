#' @title Unstabilized Inverse Probability of Treatment Weighting
#' @description Compute unstabilized IPTW for time-varying treatments.
#' @name unstabilized_iptw
#' @param identifier Name of the column for unique identifier.
#' @param treatment Treatment.
#' @param covariates List of time-varying covariates.
#' @param baseline Baseline covariates.
#' @param total_follow_up Total length of follow-up.
#' @param obsdata Observed data in wide format.
#' @return Unstabilized inverse of treatment probabilities
#' @noRd
#' @export
#' @author Awa Diop, Denis Talbot
#' @note This function requires data in a wide format.

unstabilized_iptw <- function(identifier, treatment, covariates, baseline, total_followup, obsdata) {
  # Check if observed data is in the correct format
  if (!is.data.frame(obsdata)) {
    stop("obsdatamust be a data frame in wide format")
  }

  # Validate the presence of required arguments
  required_args <- c("identifier", "treatment", "covariates", "baseline", "total_followup", "obsdata")
  missing_args <- setdiff(required_args, names(match.call()))
  if (length(missing_args) > 0) {
    stop(paste(missing_args, collapse = ", "), " not specified")
  }

  # Initialize matrix for temporary weights
  weights_temp <- matrix(numeric(0), nrow = nrow(obsdata), ncol = length(treatment))

  # Loop through each time point
  for (i in 2:length(treatment)) {
    # Generate variable names for treatments and covariates
    past_treatments <- paste0(treatment[1:(i-1)], collapse = "+")
    current_covariates <- paste0(unlist(covariates[1:i]), collapse = "+")
    baseline_covariates <- paste0(baseline, collapse = "+")

    # Formulating the denominator model
    form_denom <- as.formula(paste(paste0(treatment[i], "~", past_treatments),
                                   current_covariates, baseline_covariates, sep = "+"))

    # Fitting the model
    fit_denom <- glm(form_denom, family = binomial(link = "logit"), data = obsdata)

    # Predicting probabilities
    ps_denom <- predict(fit_denom, type = "response")

    # Calculating weights
    weights_temp[, i] <- (1 - obsdata[, treatment[i]]) / (1 - ps_denom) + obsdata[, treatment[i]] / ps_denom
  }

  # Calculate weights for the first time point
  first_time_covariates <- paste0(covariates[grep("1", covariates)], collapse = "+")
  form_denom_t1 <- as.formula(paste0(treatment[1], "~", paste0(unlist(covariates[1]),collapse = "+"), "+", paste0(baseline, collapse = "+")))
  fit_denom_t1 <- glm(form_denom_t1, family = binomial(link = "logit"), data = obsdata)
  ps_denom_t1 <- predict(fit_denom_t1, type = "response")
  weights_temp[, 1] <- (1 - obsdata[, treatment[1]]) / (1 - ps_denom_t1) + obsdata[, treatment[1]] / ps_denom_t1

  # Calculate cumulative product of weights
  weights <- t(apply(weights_temp, 1, cumprod))

  return(list(weights))
}
