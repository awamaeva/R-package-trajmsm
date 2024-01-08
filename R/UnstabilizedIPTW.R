#' @title Unstabilized Inverse Probability of Treatment Weighting
#' @description Compute unstabilized IPTW for time-varying treatments.
#' @name UnstabilizedIPTW
#' @param Identifier Name of the column for unique identifier.
#' @param Baseline Baseline covariates.
#' @param Covariates Time-varying covariates.
#' @param Treatment treatment.
#' @param TotalFollowUp Total length of follow-up.
#' @param ObservedData Observed data in wide format.
#' @return Unstabilized inverse of treatment probabilities
#' @noRd
#' @export UnstabilizedIPTW
#' @author Awa Diop, Denis Talbot
#' @note This function requires data in a wide format.

UnstabilizedIPTW <- function(Identifier, Baseline, Covariates, Treatment, TotalFollowUp, ObservedData) {
  # Check if observed data is in the correct format
  if (class(ObservedData) != "data.frame") {
    stop("ObservedData must be a data frame in wide format")
  }
  
  # Validate the presence of required arguments
  RequiredArgs <- c("Identifier", "Baseline", "Covariates", "Treatment", "TotalFollowUp", "ObservedData", "Time")
  MissingArgs <- setdiff(RequiredArgs, names(match.call()))
  if (length(MissingArgs) > 0) {
    stop(paste(MissingArgs, collapse = ", "), " not specified")
  }
  
  # Generate variable names for treatments and covariates
  PastTreatments <- paste0(Treatment[1:(i-1)], collapse = "+")
  CurrentCovariates <- paste0(unlist(Covariates[1:i]), collapse = "+")
  BaselineCovariates <- paste0(Baseline, collapse = "+")
  
  
  # Initialize matrix for temporary weights
  WeightsTemp <- matrix(numeric(0), nrow = nrow(ObservedData), ncol = length(TreatmentVars))
  
  # Loop through each time point
  for (i in 2:length(TreatmentVars)) {
    # Formulating the denominator model
    FormDenom <- as.formula(paste(paste0(Treatment[i], "~", PastTreatments), 
                                  CurrentCovariates ,   BaselineCovariates, sep = "+"))
    
    # Fitting the model
    FitDenom <- glm(FormDenom, family = binomial(link = "logit"), data = ObservedData)
    
    # Predicting probabilities
    PSDenom <- predict(FitDenom, type = "response")
    
    # Calculating weights
    WeightsTemp[, i] <- (1 - ObservedData[, Treatment[i]]) / (1 - PSDenom) + ObservedData[, Treatment[i]] / PSDenom
  }
  
  # Calculate weights for the first time point
  FirstTimeCovariates <- paste0(Covariates[grep("1", Covariates)], collapse = "+")
  FormDenomT1 <- as.formula(paste0(Treatment[1], "~", FirstTimeCovariates, "+", paste0(Baseline, collapse = "+")))
  FitDenomT1 <- glm(FormDenomT1, family = binomial(link = "logit"), data = ObservedData)
  PSDenomT1 <- predict(FitDenomT1, type = "response")
  WeightsTemp[, 1] <- (1 - ObservedData[, Treatment[1]]) / (1 - PSDenomT1) + ObservedData[, Treatment[1]] / PSDenomT1
  
  # Calculate cumulative product of weights
  Weights <- t(apply(WeightsTemp, 1, cumprod))
  
  # Return the computed weights
  return(Weights)
}

