#' @title Stabilized Inverse Probability of Treatment Weighting
#' @description Compute Stabilized IPTW for time-varying treatments.
#' @name StabilizedIPTW
#' @param Identifier Name of the column for unique identifier.
#' @param Treatment Time-varying treatment.
#' @param Covariates Time-varying covariates.
#' @param Baseline Baseline covariates.
#' @param TotalFollowUp Total length of follow-up.
#' @param ObservedData Observed data in wide format.
#' @return Stabilized inverse of treatment probabilities
#' @noRd
#' @export StabilizedIPTW
#' @author Awa Diop, Denis Talbot

StabilizedIPTW <- function(Identifier, Treatment, Covariates, Baseline, TotalFollowUp, ObservedData) {
  # Ensure observed data is a data frame
  if (class(ObservedData) != "data.frame") {
    stop("ObservedData needs to be a data frame")
  }
  
  # Check for missing required arguments
  RequiredArgs <- c("Identifier", "Baseline", "Covariates", "Treatment", "TotalFollowUp", "ObservedData")
  MissingArgs <- setdiff(RequiredArgs, names(match.call()))
  if (length(MissingArgs) > 0) {
    stop(paste(MissingArgs, collapse = ", "), " not specified")
  }
  
  # Initialize the matrix for temporary stabilized weights
  SWTemp <- matrix(numeric(0), nrow = nrow(ObservedData), ncol = length(Treatment))
  
  # Loop through each treatment time point
  for (i in 2:length(Treatment)) {
    PastTreatments <- paste0(Treatment[1:(i-1)], collapse = "+")
    CurrentCovariates <- paste0(unlist(Covariates[1:i]), collapse = "+")
    BaselineCovariates <- paste0(Baseline, collapse = "+")
    
    # Formulating the numerator and denominator models
    FormNum <- as.formula(paste0(Treatment[i], "~", PastTreatments))
    FormDenom <- as.formula(paste(Treatment[i], "~", PastTreatments, "+", CurrentCovariates, "+", BaselineCovariates))
    
    # Fitting the models
    FitNum <- glm(FormNum, family = binomial(link = "logit"), data = ObservedData)
    FitDenom <- glm(FormDenom, family = binomial(link = "logit"), data = ObservedData)
    
    # Calculating predicted probabilities
    PSNum <- ifelse(!is.na(ObservedData[, Treatment[i]]), predict(FitNum, type = "response"), NA)
    PSDenom <- ifelse(!is.na(ObservedData[, Treatment[i]]), predict(FitDenom, type = "response"), NA)
    
    # Computing the stabilized weights for each time point
    SWTemp[, i] <- ((1 - ObservedData[, Treatment[i]]) * (1 - PSNum)) / (1 - PSDenom) + (ObservedData[, Treatment[i]] * PSNum) / PSDenom
  }
  
  # Handling the first time point
  FirstTimeCovariates <- paste0(Covariates[grep("1", Covariates)], collapse = "+")
  FormNumT1 <- as.formula(paste0(Treatment[1], "~ 1"))
  FormDenomT1 <- as.formula(paste0(Treatment[1], "~", FirstTimeCovariates, "+", paste0(Baseline, collapse = "+")))
  
  FitNumT1 <- glm(FormNumT1, family = binomial(link = "logit"), data = ObservedData)
  FitDenomT1 <- glm(FormDenomT1, family = binomial(link = "logit"), data = ObservedData)
  
  PSNumT1 <- ifelse(!is.na(ObservedData[, Treatment[1]]), predict(FitNumT1, type = "response"), NA)
  PSDenomT1 <- ifelse(!is.na(ObservedData[, Treatment[1]]), predict(FitDenomT1, type = "response"), NA)
  
  # Compute the stabilized weight for the first time point
  SWTemp[, 1] <- ((1 - ObservedData[, Treatment[1]]) * (1 - PSNumT1)) / (1 - PSDenomT1) + (ObservedData[, Treatment[1]] * PSNumT1) / PSDenomT1
  
  # Calculate cumulative product of stabilized weights
  Weights <- t(apply(SWTemp, 1, cumprod))
  
  # Return the computed stabilized weights
  return(Weights)
}
