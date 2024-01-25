#' @title Generate Data Trajectories for MSM
#' @param n Number of observations to generate.
#' @param include_censor Logical, if TRUE, includes censoring.
#' @param format Character, either "long" or "wide" for the format of the output data frame.
#' @export gendata
#' @return A data frame with generated trajectories.
#' @examples
#' gendata(n = 100, include_censor = FALSE, format = "wide",total_followup = 3)

gendata<- function(n, include_censor = FALSE, format = c("long", "wide"), total_followup, seed) {
  set.seed(seed)

  # Common variables for all scenarios
  id <- 1:n
  age <- rbinom(n, 1, 0.5)
  sex <- rbinom(n, 1, 0.5)

  # Initialize variables
  statins <- hyper <- bmi <- censor <- matrix(NA, nrow = n, ncol = total_followup)
  bmi[, 1] <- rbinom(n, 1, plogis(0.5 * age + 0.7 * sex))
  hyper[, 1] <- rbinom(n, 1, plogis(.25 * age + 0.7 * sex + 0.1 * bmi[, 1]))
  statins[, 1] <- rbinom(n, 1, plogis(.75 + 0.4 * age + 0.25* sex - 0.1 * bmi[, 1] - 0.2*hyper[, 1]))
  censor[, 1] <-  rbinom(n, 1, plogis(-2 + 0.2 * age + 0.01 * sex + 0.1 * bmi[, 1] - 0.2*hyper[, 1] - 0.5*statins[, 1]))
  # Generate data based on conditions
  for (i in 2:total_followup) {
    bmi[, i] <- rbinom(n, 1, plogis(0.5 * age + 0.7 * sex - 0.65 * statins[, i-1]))
    hyper[, i] <- rbinom(n, 1, plogis(0.5 * age + 0.7 * sex + 0.1 * bmi[, i] - 0.75 * statins[, i-1]))
    statins[, i] <- rbinom(n, 1, plogis(.5 + 0.1 * age + 0.1 * sex - 0.1 * bmi[, i] - 0.2*hyper[, i] + 0.75 * statins[, i-1] ))
    if (include_censor) {
      censor[, i] <- rbinom(n, 1, plogis(-2 + 0.02 * age + 0.01 * sex - 0.5 * statins[, i] + 0.1 * hyper[, i] + 0.2 * bmi[, i]))
    }
  }

  # Apply censoring
  if (include_censor) {
    for (i in 2:total_followup) {
      statins[, i][censor[, i-1] == 1] <- NA
      hyper[, i][censor[, i-1] == 1] <- NA
      bmi[, i][censor[, i-1] == 1] <- NA
      censor[, i][censor[, i-1] == 1] <- 1
    }
  }

  y <- rbinom(n, 1, plogis(0.05 - rowSums(sapply(1:3, function(i) c(-0.25, -0.5, -0.75)[i] * statins[, i] + c(0.2, 0.3, 0.4)[i] * hyper[, i] + c(0.1, 0.2, 0.3)[i] * bmi[, i]), na.rm = TRUE)))

  obsdata <- data.frame(id = id, age = age, sex = sex,statins, hyper, bmi,y = y)

  coln <- paste0(rep(c("statins", "hyper", "bmi"), each = total_followup), seq(2011, 2010+total_followup, by = 1))
  colnames(obsdata)[4:(3+length(coln))] <- coln

  if (format == "long") {
    obsdata_long <- reshape(obsdata, varying = list(coln),
                                     v.names = c("statins", "hyper", "bmi"),
                                  idvar = "id", timevar = "time", times = rep(seq(2011, 2010+total_followup, by = 1),3), direction = "long")
    if (include_censor) {
      obsdata_long$censor <- as.vector(censor)
    }
    return(obsdata_long)
  }

  if (include_censor && format == "wide") {
    colnames(censor) <- paste0("censor", seq(2011, 2010+total_followup, by = 1))
    obsdata <- cbind(obsdata, censor)
  }

  return(obsdata)
}

obsdata <- gendata(n = 100, format = "wide", total_followup = 5, seed = 345)
