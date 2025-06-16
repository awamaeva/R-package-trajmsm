#' @title Inverse Probability of Censoring Weights
#' @description Compute unstabilized IPCW for time-varying treatments.
#' @name unstabilized_ipcw
#' @param identifier Name of the column for unique identifier.
#' @param treatment Time-varying treatment.
#' @param covariates Time-varying covariates.
#' @param baseline Baseline covariates.
#' @param censor Name of the censoring variable.
#' @param obsdata Observed data in wide format.
#' @return Unstabilized Inverse Probability of Censoring Weights
#' @keywords internal
#' @noRd
#' @author Awa Diop, Denis Talbot
#' @note This function requires data in a wide format.
#' @examples
#' obsdata = gendata(n = 1000, format = "wide",include_censor = TRUE,seed = 945)
#' baseline_var <- c("age","sex")
#' covariates <- list(c("hyper2011", "bmi2011"),
#' c("hyper2012", "bmi2012"),c("hyper2013", "bmi2013"))
#' treatment_var <- c("statins2011","statins2012","statins2013")
#' censor_var <- c("censor2011", "censor2012","censor2013")
#' unstabilized_weights = unstabilized_ipcw(
#' identifier = "id", covariates = covariates, treatment = treatment_var,
#' baseline = baseline_var, censor = censor_var,obsdata = obsdata)
#' summary(unstabilized_weights[[1]])

unstabilized_ipcw <- function(identifier, treatment, covariates, baseline, censor, obsdata) {
  if (!is.data.frame(obsdata)) {
    stop("obsdata must be a data frame in wide format")
  }

  required_args <- c("identifier", "baseline", "covariates", "treatment", "censor", "obsdata")
  missing_args <- setdiff(required_args, names(match.call()))
  if (length(missing_args) > 0) {
    stop(paste(missing_args, collapse = ", "), " not specified")
  }

  weight_treatment_temp <- matrix(numeric(0), nrow = nrow(obsdata), ncol = length(treatment))
  weight_censor_temp <- matrix(numeric(0), nrow = nrow(obsdata), ncol = length(treatment))

  for (i in 2:length(treatment)) {
    past_treatments <- paste0(treatment[1:(i-1)], collapse = "+")
    current_covariates <- paste0(unlist(covariates[1:i]), collapse = "+")
    baseline_covariates <- paste0(baseline, collapse = "+")

    form_denom_treatment <- as.formula(paste(paste0(treatment[i], "~", past_treatments, "+", current_covariates, "+", baseline_covariates)))
    form_denom_censor <- as.formula(paste(paste0(censor[i], "~", past_treatments, "+", current_covariates, "+", baseline_covariates)))

    fit_denom_treatment <- glm(form_denom_treatment, family = binomial(link = "logit"), data = obsdata)
    fit_denom_censor <- glm(form_denom_censor, family = binomial(link = "logit"), data = obsdata)

    ps_denom_treatment <- predict(fit_denom_treatment, type = "response")
    ps_denom_censor <- predict(fit_denom_censor, type = "response")

    weight_treatment_temp[, i][obsdata[,censor[i-1]] == 0]  <- (obsdata[, treatment[i]][obsdata[,censor[i-1]] == 0] / ps_denom_treatment[obsdata[,censor[i-1]] == 0] ) + ((1 - obsdata[, treatment[i]][obsdata[,censor[i-1]] == 0]) / (1 - ps_denom_treatment[obsdata[,censor[i-1]] == 0]))
    weight_censor_temp[, i][obsdata[,censor[i-1]] == 0] <- 1 / (1 - ps_denom_censor)
  }

  form_denom_treatment_t1 <- as.formula(paste0(treatment[1], "~", paste0(unlist(covariates[1]),collapse = "+"), "+", paste0(baseline, collapse = "+")))
  form_denom_censor_t1 <- as.formula(paste0(censor[1], "~", paste0(unlist(covariates[1]),collapse = "+"), "+", paste0(baseline, collapse = "+")))

  fit_denom_treatment_t1 <- glm(form_denom_treatment_t1, family = binomial(link = "logit"), data = obsdata)
  fit_denom_censor_t1 <- glm(form_denom_censor_t1, family = binomial(link = "logit"), data = obsdata)

  ps_denom_treatment_t1 <- predict(fit_denom_treatment_t1, type = "response")
  ps_denom_censor_t1 <- predict(fit_denom_censor_t1, type = "response")

  weight_treatment_temp[, 1] <- (obsdata[, treatment[1]] / ps_denom_treatment_t1) + ((1 - obsdata[, treatment[1]]) / (1 - ps_denom_treatment_t1))
  weight_censor_temp[, 1] <- 1 / (1 - ps_denom_censor_t1)

  weight_censor_treatment <- weight_treatment_temp * weight_censor_temp

  weights <- t(apply(weight_censor_treatment, 1, cumprod))

  return(list(weights))
}
