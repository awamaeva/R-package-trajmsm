#' @title Inverse Probability Weighting
#' @description Compute stabilized and unstabilized with and without censor weights.
#' @name InverseProbabilityWeighting
#' @param Identifier Name of the column for unique identifier.
#' @param Baseline Baseline covariates.
#' @param Covariates Time-varying covariates.
#' @param Treatment Time-varying treatment.
#' @param IncludeCensor Logical value TRUE/FALSE to include or not a censoring variable.
#' @param TotalFollowUp Total length of follow-up.
#' @param CensorVariable Name of the censoring variable.
#' @param ObservedData Observed data in wide format.
#' @return Inverse Probability Weights (Stabilized and Unstabilized) with and without censoring
#' @author Awa Diop, Denis Talbot
#' @export InverseProbabilityWeighting
#' @examples
#' dat = gendatTrajMSM(n = 500, Censor = FALSE)
#' Baseline <- c("Age", "Sex")
#' Covariates <- c("Hyper", "BMI")
#' Treatment <- "Statins"
#' CensorVariable <- "C"
#' sw = IPW(numerator = "stabilized", Identifier = "ID", Baseline = c("Age", "Sex"),
#'         Covariates = c("Hyper", "BMI"), Time = "Time",
#'         Treatment = "Statins", ObservedData = dat)

InverseProbabilityWeighting <- function(numerator = c("stabilized", "unstabilized"), Identifier = id,
                Baseline = V, Covariates = L, Treatment = A, IncludeCensor = FALSE,
                TotalFollowUp = K, CensorVariable = censor, 
                ObservedData = obsdata) {
  
  # Adjust function calls based on the 'numerator' and 'IncludeCensor' parameters
  if (numerator == "unstabilized" && !IncludeCensor) {
    IPWResult = UnstabilizedIPTW(Identifier, Baseline, Covariates, Treatment, TotalFollowUp, ObservedData)
  }
  
  if (numerator == "stabilized" && !IncludeCensor) {
    IPWResult = StabilizedIPTW(Identifier, Baseline, Covariates, Treatment, TotalFollowUp, ObservedData)
  }
  
  if (numerator == "unstabilized" && IncludeCensor) {
    IPWResult = UnstabilizedIPCW(Identifier, Baseline, Covariates, Treatment, TotalFollowUp, CensorVariable, ObservedData)
  }
  
  if (numerator == "stabilized" && IncludeCensor) {
    IPWResult = StabilizedIPCW(Identifier, Baseline, Covariates, Treatment, TotalFollowUp, CensorVariable, ObservedData)
  }
  
  class(IPWResult) <- "InverseProbabilityWeighting"
  return(IPWResult)
}
