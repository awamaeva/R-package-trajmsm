#' @title Generate Data Trajectories for MSM
#' @param n Number of observations to generate.
#' @param IncludeCensor Logical, if TRUE, includes censoring.
#' @param Format Character, either "long" or "wide" for the format of the output data frame.
#' @export GenDataTrajMSM
#' @return A data frame with generated trajectories.
#' @examples
#' GenDataTrajmsm(n = 100, IncludeCensor = FALSE, format = "wide")
GenDataTrajMSM <- function(n, IncludeCensor = FALSE, Format = c("long", "wide"), Seed) {

   set.seed(Seed)

  # Common variables for all scenarios
  ID = 1:n
  Time = rep(c(2011, 2012, 2013), n)
  Age = rbinom(n, 1, 0.5)
  Sex = rbinom(n, 1, 0.5)

  # Initialize variables
  Statins <- Hyper <- BMI <- Censor <- matrix(NA, nrow = n, ncol = 3)

  # Generate data based on conditions
  Statins[, 1] = rbinom(n, 1, plogis(0.25 + 0.2 * Age + 0.1 * Sex))
  BMI[, 1] = rbinom(n, 1, plogis(0.5 * Age + 0.7 * Sex - 0.5 * Statins[, 1]))
  Hyper[, 1] = rbinom(n, 1, plogis(0.5 * Age + 0.7 * Sex + 0.1 * BMI[, 1]- 0.5 * Statins[, 1]))
  Censor[, 1] = rbinom(n, 1, plogis(-2.5 - 0.5 * Statins[, 1] + 0.5 * Hyper[, 1] + 0.2 * BMI[, 1]))

  for (i in 2:3) {
    Hyper[, i] = rbinom(n, 1, plogis(0.5 * Age + 0.7 * Sex - 0.5 * Statins[, i-1] + Hyper[, i-1] + 0.1 * BMI[, i-1]))
    BMI[, i] = rbinom(n, 1, plogis(0.5 * Age + 0.7 * Sex - 0.1 * Statins[, i-1] +  0.2 * Hyper[, i-1] + BMI[, i-1]))
    Statins[, i] = rbinom(n, 1, plogis(-0.1 + 0.2 * Age + 0.1 * Sex + 0.5 * Statins[, i-1] + 0.2 * Hyper[, i] +  0.1 * BMI[, i]))

    if (IncludeCensor) {
       Censor[, i] = rbinom(n, 1, plogis(-2.5 - 0.5 * Statins[, i-1] + 0.5 * Hyper[, i-1] + 0.2 * BMI[, i-1]))
        Censor[, i][Censor[, i-1] == 1] <- 1
    }

}

  for (i in 2:3) {
    if (IncludeCensor) {
    Statins[, i][Censor[, i-1] == 1] <- NA
    Hyper[, i][Censor[, i-1] == 1] <- NA
    BMI[, i][Censor[, i-1] == 1] <- NA
    }
  }

  Y = rbinom(n, 1, plogis(0.05 - rowSums(sapply(1:3, function(i) c(-0.1, -0.15, -0.2)[i] * Statins[, i] + c(0.2, 0.3, 0.4)[i] * Hyper[, i] + c(0.1, 0.2, 0.3)[i] * BMI[, i]), na.rm = TRUE)))

  ObservedData <- data.frame(ID = ID, Age = Age, Sex = Sex,
                             Statins2011 = Statins[, 1], Statins2012 = Statins[, 2], Statins2013 = Statins[, 3],
                             Hyper2011 = Hyper[, 1], Hyper2012 = Hyper[, 2], Hyper2013 = Hyper[, 3],
                             BMI2011 = BMI[, 1], BMI2012 = BMI[, 2], BMI2013 = BMI[, 3], Y = Y)
  # Reshape data based on format
  if (Format == "long") {
    ObservedData <- data.frame(ID = rep(ID, 3), Time = Time, Age = rep(Age, 3), Sex = rep(Sex, 3),
                      Statins = as.vector(Statins), Hyper = as.vector(Hyper), BMI = as.vector(BMI), Y = rep(Y, 3))}

  if (IncludeCensor & Format == "long"){
    ObservedData$Censor = as.vector(Censor)
  }

    if (IncludeCensor & Format == "wide"){

      ObservedData$Censor2011 = Censor[, 1]
      ObservedData$Censor2012 = Censor[, 2]
      ObservedData$Censor2013 = Censor[, 3]
    }

  return(ObservedData)
}

ObsData <- GenDataTrajMSM(n = 100, IncludeCensor = TRUE, Format = "wide", Seed = 345)
