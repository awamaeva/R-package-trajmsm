#' @title Stabilized Inverse Probability of Treatment Weighting
#' @description Compute Stabilized IPTW for time-varying treatments.
#' @name stabilized_iptw
#' @param identifier Name of the column for unique identifier.
#' @param treatment Time-varying treatment.
#' @param covariates Time-varying covariates.
#' @param baseline Baseline covariates.
#' @param obsdata Observed data in wide format.
#' @return Stabilized inverse of treatment probabilities
#' @keywords internal
#' @noRd
#' @note This function requires data in a wide format.
#' @author Awa Diop, Denis Talbot
#' @examples
#' obsdata = gendata(n = 1000, format = "wide",total_followup =3, seed = 945)
#' baseline_var <- c("age","sex")
#' covariates <- list(c("hyper2011", "bmi2011"),
#' c("hyper2012", "bmi2012"),c("hyper2013", "bmi2013"))
#' treatment_var <- c("statins2011","statins2012","statins2013")
#' stabilized_weights = stabilized_iptw(
#' identifier = "id", covariates = covariates, treatment = treatment_var,
#' baseline = baseline_var,obsdata = obsdata)
#' summary(stabilized_weights[[1]])


stabilized_iptw <- function(identifier, treatment, covariates, baseline, obsdata) {
  if (!is.data.frame(obsdata)) {
    stop("obsdata needs to be a data frame")
  }

  required_args <- c("identifier", "baseline", "covariates", "treatment", "obsdata")
  missing_args <- setdiff(required_args, names(match.call()))
  if (length(missing_args) > 0) {
    stop(paste(missing_args, collapse = ", "), " not specified")
  }

  sw_temp <- matrix(numeric(0), nrow = nrow(obsdata), ncol = length(treatment))

  for (i in 2:length(treatment)) {
    past_treatments <- paste0(treatment[1:(i-1)], collapse = "+")
    current_covariates <- paste0(unlist(covariates[1:i]), collapse = "+")
    baseline_covariates <- paste0(baseline, collapse = "+")

    form_num <- as.formula(paste0(treatment[i], "~", past_treatments))
    form_denom <- as.formula(paste(treatment[i], "~", past_treatments, "+", current_covariates, "+", baseline_covariates))

    fit_num <- glm(form_num, family = binomial(link = "logit"), data = obsdata)
    fit_denom <- glm(form_denom, family = binomial(link = "logit"), data = obsdata)

    ps_num <- ifelse(!is.na(obsdata[, treatment[i]]), predict(fit_num, type = "response"), NA)
    ps_denom <- ifelse(!is.na(obsdata[, treatment[i]]), predict(fit_denom, type = "response"), NA)

    sw_temp[, i] <- ((1 - obsdata[, treatment[i]]) * (1 - ps_num)) / (1 - ps_denom) + (obsdata[, treatment[i]] * ps_num) / ps_denom
  }

  first_time_covariates <- paste0(covariates[grep("1", covariates)], collapse = "+")
  form_num_t1 <- as.formula(paste0(treatment[1], "~ 1"))

  form_denom_t1 <- as.formula(paste0(treatment[1], "~", paste0(unlist(covariates[1]),collapse = "+"), "+", paste0(baseline, collapse = "+")))

  fit_num_t1 <- glm(form_num_t1, family = binomial(link = "logit"), data = obsdata)
  fit_denom_t1 <- glm(form_denom_t1, family = binomial(link = "logit"), data = obsdata)

  ps_num_t1 <- ifelse(!is.na(obsdata[, treatment[1]]), predict(fit_num_t1, type = "response"), NA)
  ps_denom_t1 <- ifelse(!is.na(obsdata[, treatment[1]]), predict(fit_denom_t1, type = "response"), NA)

  sw_temp[, 1] <- ((1 - obsdata[, treatment[1]]) * (1 - ps_num_t1)) / (1 - ps_denom_t1) + (obsdata[, treatment[1]] * ps_num_t1) / ps_denom_t1

  weights <- t(apply(sw_temp, 1, cumprod))

  return(list(weights))
}
