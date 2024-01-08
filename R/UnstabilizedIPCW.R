#' @title Inverse Probability of Censoring Weights
#' @description Compute unstabilized IPCW for time-varying treatments.
#' @name UnstabilizedIPCW
#' @param Identifier Name of the column for unique identifier.
#' @param Baseline Baseline covariates.
#' @param Covariates Time-varying covariates.
#' @param Treatment Time-varying treatment.
#' @param TotalFollowUp Total length of follow-up.
#' @param CensorVariable Name of the censoring variable.
#' @param Time Time variable.
#' @param ObservedData Observed data in wide format.
#' @return Unstabilized Inverse Probability of Censoring Weights
#' @noRd
#' @export UnstabilizedIPCW
#' @author Awa Diop, Denis Talbot
#' @note This function requires data in a wide format.

UnstabilizedIPCW <- function(Identifier, Baseline, Covariates, Treatment, TotalFollowUp, CensorVariable, ObservedData) {
  # Validate the data format
  if (class(ObservedData) != "data.frame") {
    stop("ObservedData must be a data frame in wide format")
  }
  
  # Check for the presence of required arguments
  RequiredArgs <- c("Identifier", "Baseline", "Covariates", "Treatment", "TotalFollowUp", "CensorVariable", "Time", "ObservedData")
  MissingArgs <- setdiff(RequiredArgs, names(match.call()))
  if (length(MissingArgs) > 0) {
    stop(paste(MissingArgs, collapse = ", "), " not specified")
  }
  
  # Generate variable names for treatments, covariates, and censoring
  PastTreatments <- paste0(Treatment[1:(i-1)], collapse = "+")
  CurrentCovariates <- paste0(unlist(Covariates[1:i]), collapse = "+")
  BaselineCovariates <- paste0(Baseline, collapse = "+")
  
  # Initialize matrices for temporary weights
  WeightTreatmentTemp <- matrix(numeric(0), nrow = nrow(ObservedData), ncol = length(TreatmentVars))
  WeightCensorTemp <- matrix(numeric(0), nrow = nrow(ObservedData), ncol = length(TreatmentVars))
  
  # Loop through each time point
  for (i in 2:length(TreatmentVars)) {
    # Formulating the denominator model for treatment and censoring
    FormDenomTreatment <- as.formula(paste(paste0(TreatmentVars[i], "~", 
                                                  PastTreatments), 
                                                  CurrentCovariates, 
                                                  BaselineCovariates, sep = "+"))
    
    FormDenomCensor <- as.formula(paste(paste0(CensorVars[i], "~", 
                                               PastTreatments), 
                                               CurrentCovariates, 
                                               BaselineCovariates, sep = "+"))
    
    # Fitting models for denominator
    FitDenomTreatment <- glm(FormDenomTreatment, family = binomial(link = "logit"), data = ObservedData)
    FitDenomCensor <- glm(FormDenomCensor, family = binomial(link = "logit"), data = ObservedData)
    
    # Predicting probabilities
    PSDenomTreatment <- predict(FitDenomTreatment, type = "response")
    PSDenomCensor <- predict(FitDenomCensor, type = "response")
    
    # Calculating weights
    WeightTreatmentTemp[, i] <- (ObservedData[, Treatment[i]] / PSDenomTreatment) + ((1 - ObservedData[, Treatment[i]]) / (1 - PSDenomTreatment))
    WeightCensorTemp[, i] <- 1 / (1 - PSDenomCensor)
  }
  
  # Calculate weights for the first time point
  # Formulating the denominator model for the first time point for treatment and censoring
  FormDenomTreatmentT1 <- as.formula(paste0(Treatment[1], "~", paste0(Baseline, collapse = "+")))
  FormDenomCensorT1 <- as.formula(paste0(Censor[1], "~", paste0(Baseline, collapse = "+")))
  
  # Fitting models for denominator at t=1
  FitDenomTreatmentT1 <- glm(FormDenomTreatmentT1, family = binomial(link = "logit"), data = ObservedData)
  FitDenomCensorT1 <- glm(FormDenomCensorT1, family = binomial(link = "logit"), data = ObservedData)
  
  # Predicting probabilities at t=1
  PSDenomTreatmentT1 <- predict(FitDenomTreatmentT1, type = "response")
  PSDenomCensorT1 <- predict(FitDenomCensorT1, type = "response")
  
  # Calculating weights for the first time point
  WeightTreatmentTemp[, 1] <- (ObservedData[, TreatmentVars[1]] / PSDenomTreatmentT1) + 
                              ((1 - ObservedData[, TreatmentVars[1]]) / (1 - PSDenomTreatmentT1))
  WeightCensorTemp[, 1] <- 1 / (1 - PSDenomCensorT1)
  
  # Compute final unstabilized weights
  WeightCensorTreatment <- WeightTreatmentTemp * WeightCensorTemp
  
  # Calculate cumulative product of weights
  Weights <- t(apply( WeightCensorTreatment, 1, cumprod))
  
  
  # Return the computed weights
  return(Weights)
}
  