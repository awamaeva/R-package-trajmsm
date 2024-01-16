#' @title Counterfactual means for g-formula.
#' @description Get the counterfactual means for the g-formula.
#' @name gformula
#' @param formula specification of the model for the outcome to be fitted.
#' @param baseline name of the baseline covariates.
#' @param covariates list of the names of the time-varying covariates.
#' @param treatment name of the time-varying treatment.
#' @param Y outcome variable.
#' @param ntimes_interval length of a time-interval (s).
#' @param obsdata observed data in wide format.
#' @returns  \item{list_gform_countermeans}{Counterfactual means obtained with g-formula.}
#' @import e1071
#' @author Awa Diop, Denis Talbot
#' @examples
#' Obsdata = gendata_trajmsm(n = 1000, format = "wide", seed = 745)
#' baseline_var <- c("age","sex")
#' covariates <- list(c("hyper2011", "bmi2011"),c("hyper2012", "bmi2012"),c("hyper2013", "bmi2013"))
#' treatment_var <- c("statins2011","statins2012","statins2013")
#' res_gform = gformula(formula = as.formula(" y ~ statins2011 + statins2012 + statins2013 + hyper2011 + bmi2011 + hyper2012 + bmi2012 +
#'                                               hyper2013 + bmi2013 + age + sex " ), baseline = baseline, covariates = covariates,
#' treatment = treatment, outcome = outcome, ntimes_interval = 3, obsdata = obsdata)
#' res_gform$counter_means

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

