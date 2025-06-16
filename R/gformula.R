#' @title Counterfactual means via G-Formula
#' @description Calculates counterfactual means using the g-formula approach.
#' @name gformula
#' @param formula Specification of the model for the outcome to be fitted.
#' @param baseline Names of the baseline covariates.
#' @param covariates Names of the time-varying covariates (should be a list).
#' @param treatment Names of the time-varying treatment.
#' @param outcome Name of the outcome variable.
#' @param ntimes_interval Length of a time-interval (s).
#' @param obsdata Observed data in wide format.
#' @returns \item{list_gform_countermeans}{ List of counterfactual means obtained with g-formula.}
#' @import e1071
#' @export
#' @examples
#' \donttest{
#' obsdata = gendata(n = 1000, format = "wide", total_followup = 6, seed = 945)
#' years <- 2011:2016
#' baseline_var <- c("age","sex")
#' variables <- c("hyper", "bmi")
#' var_cov <- c("statins","hyper", "bmi")
#' covariates <- lapply(years, function(year) {
#' paste0(variables, year)})
#' treatment_var <- paste0("statins", 2011:2016)
#' formula = paste0("y ~", paste0(treatment_var,collapse = "+"), "+",
#'                 paste0(unlist(covariates), collapse = "+"),"+",
#'                 paste0(baseline_var, collapse = "+"))
#'res_gform <- gformula(formula = formula, baseline = baseline_var, covariates = covariates,
#'treatment = treatment_var, outcome = "y", ntimes_interval = 6, obsdata =   obsdata )
#'}

#' @author Awa Diop, Denis Talbot


gformula <- function(formula, baseline, covariates, treatment, outcome, ntimes_interval, obsdata) {
nregimes <- 2^ntimes_interval  # Number of treatment regimes

  # Treatment regimes
  dat_combinations <- bincombinations(ntimes_interval)
  list_regimes <- lapply(1:nregimes, function(x) {
    regime <- dat_combinations[x, ]
    return(regime)
  })

  # Matrix to store counterfactual means
  res_counter_means <- matrix(0, nrow = nregimes, ncol = 1)
  j <- 0  # Index for storing results
  mod_Qsp1 <- glm(formula = formula, family = binomial, data = obsdata)

  # For each treatment regime
  for (regime in list_regimes) {
    j <- j + 1  # Update regime index
    dat0 <- obsdata
    dat0[, treatment[1:ntimes_interval]] <- do.call(cbind, lapply(1:ntimes_interval, function(x) {
      rgm <- rep(regime[x], nrow(dat0))
      return(rgm)
    }))

    obsdata$Qs<- predict(mod_Qsp1, newdata = dat0, type = "response")

    for (i in (ntimes_interval - 1):1) {

      # Update formula
      current_treatments <- paste0(treatment[1:i], collapse = "+")
      current_covariates <- paste0(unlist(covariates[1:i]), collapse = "+")
      baseline_covariates <- paste0(baseline, collapse = "+")

      formula_updated <- as.formula(paste0("Qs~ ", current_treatments, " + ", current_covariates, " + ", baseline_covariates))

      # Estimation
      mod_Qs <- glm(formula_updated, family = binomial, data = obsdata)

      # Prediction under a fixed treatment regime
      dat0 <- obsdata
      for (k in 1:i) {
        treatment_variable <- treatment[k]
        dat0[[treatment_variable]] <- as.numeric(rep(regime[k], nrow(dat0)))
      }

      obsdata$Qs <- predict(mod_Qs, newdata = dat0, type = "response")
    }
    Q0 <- mean(obsdata$Qs, na.rm = TRUE)

    res_counter_means[j, ] <- Q0

  }

  list_gform_countermeans <- list(treatment_regimes = list_regimes, counter_means = res_counter_means)

  return(list_gform_countermeans)
}

