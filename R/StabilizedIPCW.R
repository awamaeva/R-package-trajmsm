#' @title Inverse Probability of Censoring Weights
#' @description Compute stabilized Inverse Probability of Censoring Weights
#' @name StabilizedIPCW
#' @param Identifier Name of the column for unique identifier.
#' @param Treatment Time-varying treatment.
#' @param Covariates Time-varying covariates.
#' @param Baseline Baseline covariates.
#' @param TotalFollowUp Total length of follow-up.
#' @param CensorVariable Name of the censoring variable.
#' @param ObservedData Observed data in wide format.
#' @return Stabilized Inverse Probability of Censoring Weights
#' @noRd
#' @export StabilizedIPCW
#' @author Awa Diop, Denis Talbot
#' @note This function requires data in a wide format.

StabilizedIPCW <- function(Identifier, Treatment, Covariates, Baseline, TotalFollowUp, CensorVariable, ObservedData) {
  # Check if observed data is in the correct format
  if (class(ObservedData) != "data.frame") {
    stop("ObservedData must be a data frame in wide format")
  }
  
  # Validate the presence of required arguments
  RequiredArgs <- c("Identifier", "Baseline", "Covariates", "Treatment", "TotalFollowUp", "CensorVariable", "ObservedData")
  MissingArgs <- setdiff(RequiredArgs, names(match.call()))
  if (length(MissingArgs) > 0) {
    stop(paste(MissingArgs, collapse = ", "), " not specified")
  }
  
  # Initialize matrices for temporary weights
  SWTreatmentTemp <- matrix(numeric(0), nrow = nrow(ObservedData), ncol = length(Treatment))
  SWCensorTemp <- matrix(numeric(0), nrow = nrow(ObservedData), ncol = length(Treatment))
  
  # Loop through each time point
  for (i in 2:length(Treatment)) {
    PastTreatments <- paste0(Treatment[1:(i-1)], collapse = "+")
    CurrentCovariates <- paste0(unlist(Covariates[1:i]), collapse = "+")
    BaselineCovariates <- paste0(Baseline, collapse = "+")
    
    # Formulating numerator model for treatment and censoring
    FormNumTreatment <- as.formula(paste(paste0(Treatment[i], "~", PastTreatments)))
    FormNumCensor <- as.formula(paste(paste0(CensorVariable[i], "~ 1")))
    
    # Formulating denominator model for treatment and censoring
    FormDenomTreatment <- as.formula(paste(paste0(Treatment[i], "~", PastTreatments, "+", 
                                                  CurrentCovariates, "+", BaselineCovariates)))
    FormDenomCensor <- as.formula(paste(paste0(CensorVariable[i], "~", PastTreatments, "+", 
                                                 CurrentCovariates, "+",  BaselineCovariates)))
    
    # Fitting models for numerator
    FitNumTreatment <- glm(FormNumTreatment, family = binomial(link = "logit"), data = ObservedData)
    FitNumCensor <- glm(FormNumCensor, family = binomial(link = "logit"), data = ObservedData)
    
    # Fitting models for denominator
    FitDenomTreatment <- glm(FormDenomTreatment, family = binomial(link = "logit"), data = ObservedData)
    FitDenomCensor <- glm(FormDenomCensor, family = binomial(link = "logit"), data = ObservedData)
    
    # Predicting probabilities for treatment and censoring
    PSTreatment <- predict(FitNumTreatment, type = "response")
    PSCensor <- predict(FitNumCensor, type = "response")
    
    # Calculating stabilized weights
    SWTreatmentTemp[, i] <- (ObservedData[, Treatment[i]] * PSTreatment) / predict(FitDenomTreatment, type = "response")
    SWCensorTemp[, i] <- (1 - PSCensor) / predict(FitDenomCensor, type = "response")
  }
  
  
  # Initialize matrices for temporary weights
  SWTreatmentTemp <- matrix(numeric(0), nrow = nrow(ObservedData), ncol = length(Treatment))
  SWCensorTemp <- matrix(numeric(0), nrow = nrow(ObservedData), ncol = length(Treatment))
  
  # Calculate stabilized weights for the first time point
  # Formulating numerator model for treatment and censoring at t=1
  FormNumTreatmentT1 <- as.formula(paste0(Treatment[1], "~ 1"))
  FormNumCensorT1 <- as.formula(paste0(CensorVariable[1], "~ 1"))
  
  # Formulating denominator model for treatment and censoring at t=1
  # Including only the baseline covariates since it's the first time point
  FormDenomTreatmentT1 <- as.formula(paste0(Treatment[1], "~", paste0(Baseline, collapse = "+")))
  FormDenomCensorT1 <- as.formula(paste0(CensorVariable[1], "~", paste0(Baseline, collapse = "+")))
  
  # Fitting models for numerator at t=1
  FitNumTreatmentT1 <- glm(FormNumTreatmentT1, family = binomial(link = "logit"), data = ObservedData)
  FitNumCensorT1 <- glm(FormNumCensorT1, family = binomial(link = "logit"), data = ObservedData)
  
  # Fitting models for denominator at t=1
  FitDenomTreatmentT1 <- glm(FormDenomTreatmentT1, family = binomial(link = "logit"), data = ObservedData)
  FitDenomCensorT1 <- glm(FormDenomCensorT1, family = binomial(link = "logit"), data = ObservedData)
  
  # Predicting probabilities for treatment and censoring at t=1
  PSTreatmentT1 <- predict(FitNumTreatmentT1, type = "response")
  PSCensorT1 <- predict(FitNumCensorT1, type = "response")
  
  # Calculating stabilized weights for the first time point
  SWTreatmentTemp[, 1] <- (ObservedData[, Treatment[1]] * PSTreatmentT1) / predict(FitDenomTreatmentT1, type = "response")
  SWCensorTemp[, 1] <- (1 - PSCensorT1) / predict(FitDenomCensorT1, type = "response")
  
  # Compute final stabilized weights
  SWCensorTreatment <- matrix(numeric(0), nrow = nrow(ObservedData), ncol = length(Treatment))
  
  # Multiply the temporary treatment and censoring weights to get the combined weights
  SWCensorTreatment <- SWTreatmentTemp * SWCensorTemp
  
  # Calculate cumulative product of weights
  Weights <- t(apply(SWCensorTreatment, 1, cumprod))
  
  # Return the computed cumulative weights
  return(Weights)
}
