#' @title Stabilized Inverse of Probability of Treatment Weighting
#' @description Compute Stabilized IPTW for time-varying treatments.
#' @name IPTW_sw
#' @param id  name of the column for unique identifiant.
#' @param V baseline covariates.
#' @param L time-varying covariates.
#' @param A time-varying treatment.
#' @param K number of measuring times.
#' @param obsdata observed data in wide format.
#' @return \item{SW}{Stabilized inverse of treatment probabilities}
#' @noRd
#' @export IPTW_sw
#' @author Awa Diop, Denis Talbot

IPTW_sw <- function(id, A, L, V, K, time, obsdata) {
  if (class(obsdata) != "data.frame") stop("Convert obsdata into a data.frame")
  
  required_args <- c("id", "V", "L", "A", "K", "time", "obsdata")
  missing_args <- setdiff(required_args, names(match.call()))
  if (length(missing_args) > 0) stop(paste(missing_args, collapse=", "), " not specified")
  
  SW_temp <- matrix(numeric(0), nrow = nrow(obsdata), ncol = length(A))
  
  for (i in 2:length(A)) {
    past_treatments <- paste0(A[1:(i-1)], collapse = "+")
    current_covariates <- paste0(unlist(L[1:i]), collapse = "+")
    baseline_covariates <- paste0(V, collapse = "+")
    
    form_num <- as.formula(paste0(A[i], "~", past_treatments))
    form_denom <- as.formula(paste(A[i], "~", past_treatments, "+", current_covariates, "+", baseline_covariates))
    
    fit_num <- glm(form_num, family = binomial(link = "logit"), data = obsdata)
    fit_denom <- glm(form_denom, family = binomial(link = "logit"), data = obsdata)
    
    ps_num <- ifelse(!is.na(obsdata[, A[i]]), predict(fit_num, type = "response"), NA)
    ps_denom <- ifelse(!is.na(obsdata[, A[i]]), predict(fit_denom, type = "response"), NA)
    
    SW_temp[, i] <- ((1 - obsdata[, A[i]]) * (1 - ps_num)) / (1 - ps_denom) + (obsdata[, A[i]] * ps_num) / ps_denom
  }
  
  # Handling the first time point
  first_time_covariates <- paste0(L[grep("1", L)], collapse = "+")
  form_num_t1 <- as.formula(paste0(A[1], "~ 1"))
  form_denom_t1 <- as.formula(paste0(A[1], "~", first_time_covariates, "+", paste0(V, collapse = "+")))
  
  fit_num_t1 <- glm(form_num_t1, family = binomial(link = "logit"), data = obsdata)
  fit_denom_t1 <- glm(form_denom_t1, family = binomial(link = "logit"), data = obsdata)
  
  ps_num_t1 <- ifelse(!is.na(obsdata[, A[1]]), predict(fit_num_t1, type = "response"), NA)
  ps_denom_t1 <- ifelse(!is.na(obsdata[, A[1]]), predict(fit_denom_t1, type = "response"), NA)
  
  SW_temp[, 1] <- ((1 - obsdata[, A[1]]) * (1 - ps_num_t1)) / (1 - ps_denom_t1) + (obsdata[, A[1]] * ps_num_t1) / ps_denom_t1
  
  SW <- t(apply(SW_temp, 1, cumprod))
  return(SW)
}
