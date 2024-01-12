#' @title Marginal Structural Model and Latent Class of Growth Analysis estimated with IPW
#' @description Combine Marginal Structural Model and Latent Class of Growth Analysis
#' @name TrajMSMIPW
#' @param Formula1 Specification of the model for the outcome to be fitted for a binomial or gaussian distribution.
#' @param Formula2 Specification of the model for the outcome to be fitted for a survival outcome.
#' @param Family Specification of the error distribution and link function to be used in the model.
#' @param Identifier Name of the column for unique identification.
#' @param Treatment Time-varying treatment.
#' @param Baseline Baseline covariates.
#' @param Covariates Time-varying covariates.
#' @param Y Outcome of interest.
#' @param TotalFollowup Number of measuring times.
#' @param NumberTraj An integer to fix the number of trajectory groups.
#' @param ObsData Dataset to be used in the analysis.
#' @param Numerator Type of weighting ("stabilized" or "unstabilized").
#' @param Weights A vector of estimated weights. If NULL, the weights are computed by the function \code{IPW}.
#' @return \item{IPW}{Stabilized and unstabilized inverse of probabilities}
#' @export
#' @importFrom sandwich
#' @importFrom survival coxph
#' @import flexmix
#' @examples
#' SampleData = GenDataTrajMSM(n = 1000, IncludeCensor = FALSE, Format = "wide", Seed = 345)
#' BaselineVar <- c("Age","Sex")
#' CoVar <- c("Hyper2011", "Hyper2012", "Hyper2013", "BMI2011", "BMI2012", "BMI2013")
#' TreatmentVar <- c("Statins2011","Statins2012","Statins2013")
#' K = 3
#' CensorVar = c("Censor2011","Censor2012","Censor2013")
#' StabilizedWeights = InverseProbabilityWeighting(Numerator = "stabilized", Identifier = "ID",
#'         Covariates = CoVar, Treatment = TreatmentVar,Baseline = BaselineVar,
#'  TotalFollowUp = K, ObsData = SampleData)
#'  LongObsData = WideToLong(ObsData = SampleData, Varying = c("Statins2011", "Statins2012", "Statins2013",
#'  "Hyper2011", "Hyper2012", "Hyper2013", "BMI2011", "BMI2012", "BMI2013"),
#'  idvar = "ID", timevar = "Time")
#' Formula = as.formula(cbind(Statins, 1 - Statins) ~ Time)
#' ResTraj = BuildTraj(ObsData = LongObsData, NumberTraj = 3, Formula = Formula, Identifier = "ID")
#' DataPost = ResTraj$DataPost
#' TrajMSMDataLong <- merge(LongObsData, DataPost, by = "ID")
#'     AggFormula <- as.formula(paste("Statins", "~", "Time", "+", "class"))
#'     AggTrajData <- aggregate(AggFormula, data = TrajMSMDataLong, FUN = mean)
#'     AggTrajData
#' TrajMSMDataWide <- merge(SampleData, DataPost, by = "ID")
#'TrajMSMDataWide[ , "TrajGroup"] <- factor(ifelse(TrajMSMDataWide[ , "class"] == "3" ,"Group1" ,
#' ifelse (TrajMSMDataWide[ , "class"]== "2" , "Group2" ,"Group3")))
#' TrajMSMDataWide[ , "TrajGroup"] <- relevel(TrajMSMDataWide[ , "TrajGroup"], ref = "Group3")
#'TrajMSMIPW(Formula1 = as.formula("Y ~ TrajGroup"),
#'            Identifier = "ID", Baseline = BaselineVar, Covariates = CoVar, Treatment = TreatmentVar,
#'            Y=Y, NumberTraj=3,TotalFollowup = 3, Family = "binomial",
#'            ObsData = TrajMSMDataWide,Numerator = "stabilized")


TrajMSMIPW <- function(Formula1, Formula2, Family, Identifier, Treatment, Covariates,
                       Baseline, Y, TotalFollowup, NumberTraj, ObsData, Numerator = "stabilized", Weights = NULL) {
  # Validate inputs
  if (!is.data.frame(ObsData)) stop("ObsData must be a data frame.")
  if (!is.character(Identifier) || !Identifier %in% names(ObsData)) stop("Identifier must be a valid column name in ObsData.")
  if (!is.character(Family) || !Family %in% c("binomial", "gaussian", "survival")) stop("Family must be 'binomial', 'gaussian', or 'survival'.")
  if (!is.character(Numerator) || !Numerator %in% c("stabilized", "unstabilized")) stop("Numerator must be 'stabilized' or 'unstabilized'.")
  if (!is.null(Weights) && length(Weights) != nrow(ObsData)) stop("Length of Weights must match the number of rows in ObsData.")

  # Compute Weights if not provided
  if (is.null(Weights)) {
    Weights <- InverseProbabilityWeighting(Identifier = Identifier, Covariates = Covariates,
                                           Treatment = Treatment, Baseline = Baseline,
                                           TotalFollowUp = TotalFollowup, Numerator = Numerator,
                                           ObsData = ObsData)[, TotalFollowup]
    ObsData$Weights <- Weights
  } else {
    ObsData$Weights <- Weights
  }

  # Model fitting
  cluster_formula <- as.formula(paste0("~", Identifier))
  if (Family == "gaussian") {
    ModGLM <- glm(Formula1, data = ObsData, weights = Weights, family = gaussian)
  } else if (Family == "binomial") {
    ModGLM <- glm(Formula1, data = ObsData, weights = Weights, family = binomial)
  } else if (Family == "survival") {
    ModGLM <- coxph(Formula2, data = ObsData, cluster = cluster_formula, weights = Weights)
  }

  # Extracting model results
  coefs <- coef(summary(ModGLM))[, 1]
  se <- sqrt(diag(vcovCL(ModGLM, cluster = cluster_formula)))
  pvalue <- coef(summary(ModGLM))[, "Pr(>|z|)"]
  IClo <- coefs - 1.96 * se
  ICup <- coefs + 1.96 * se

  ResTrajMSM <- cbind(coefs, se, pvalue, IClo, ICup)
  colnames(ResTrajMSM) <- c("Estimate", "Std.Error", "Pvalue", "Lower CI", "Upper CI")

  return(ResTrajMSM)
}
