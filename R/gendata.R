#' @title Generate data trajectories for MSM
#' @description
#' Provides datasets for running examples for LCGA-MSM and LCGA-HRMSM.
#' @name  gendata
#' @param n Number of observations to generate.
#' @param include_censor Logical, if TRUE, includes censoring.
#' @param format Character, either "long" or "wide" for the format of the output data frame.
#' @param timedep_outcome Logical, if TRUE, includes a time-dependent outcome.
#' @param start_year Baseline year.
#' @param total_followup Number of measuring times.
#' @param seed, Use a specific seed value to ensure the simulated data is replicable.
#' @importFrom stats na.omit rbinom plogis qlogis  reshape glm
#' binomial coef as.formula
#' @export
#' @return A data frame with generated data trajectories.
#' @examples
#' gendata(n = 100, include_censor = FALSE, format = "wide",total_followup = 3, seed = 945)

gendata<- function(n, include_censor = FALSE, format = c("long", "wide"),
                   start_year = 2011, total_followup, timedep_outcome = FALSE, seed) {
  set.seed(seed)

  # Common variables for all scenarios
  id <- 1:n
  age <- rbinom(n, 1, 0.5)
  sex <- rbinom(n, 1, 0.5)

  # Initialize variables
  statins <- hyper <- bmi <- censor <- matrix(NA, nrow = n, ncol = total_followup)
  bmi[, 1] <-   rbinom(n, 1, plogis(0.15 * age + 0.7 * sex))
  hyper[, 1] <-  rbinom(n, 1, plogis(.15 * age + 0.7 * sex + 0.1 * bmi[, 1]))
  statins[, 1] <- rbinom(n, 1, plogis(.5 + 0.4 * age + 0.25* sex - 0.1 * bmi[, 1] - 0.2*hyper[, 1]))
  censor[, 1] <-  rbinom(n, 1, plogis(-2 + 0.2 * age + 0.01 * sex + 0.1 * bmi[, 1] - 0.2*hyper[, 1] - 0.5*statins[, 1]))
  # Generate data based on conditions
  for (i in 2:total_followup) {
    bmi[, i] <- rbinom(n, 1, plogis(0.15 * age + 0.7 * sex - 0.25 * statins[, i-1]))
    hyper[, i] <- rbinom(n, 1, plogis(0.15 * age + 0.7 * sex + 0.1 * bmi[, i] - 0.35 * statins[, i-1]))
    statins[, i] <- rbinom(n, 1, plogis(.5 + 0.1 * age + 0.1 * sex - 0.1 * bmi[, i] - 0.2*hyper[, i] + (seq(0.15,1,length.out = total_followup)[i-1])* statins[, i-1] ))
    if (include_censor) {
      censor[, i] <- rbinom(n, 1, plogis(-2 + 0.02 * age + 0.01 * sex - 0.5 * statins[, i] + 0.1 * hyper[, i] + 0.2 * bmi[, i]))
    }
  }

  if(timedep_outcome == FALSE){
  y <- rbinom(n, 1, plogis(-2.5 + rowSums(sapply(1:total_followup, function(i)
 -0.5* statins[, i] + 0.25 * hyper[, i] + 0.25* bmi[, i]))))

  obsdata <- data.frame(id = id, age = age, sex = sex,y = y,statins = statins, hyper = hyper, bmi = bmi)

  #coln <- paste0(rep(c("statins", "hyper", "bmi"), each = total_followup), seq(start_year, start_year+total_followup-1, by = 1))
  coln <- paste0(rep(c("statins", "hyper", "bmi"), each = total_followup),
                 rep(seq(start_year, start_year + total_followup - 1), times = 3))

  colnames(obsdata)[5:ncol(obsdata)] <- coln

  # Apply censoring
  if (include_censor) {
    for (i in 2:total_followup) {
      statins[, i][censor[, i-1] == 1] <- NA
      hyper[, i][censor[, i-1] == 1] <- NA
      bmi[, i][censor[, i-1] == 1] <- NA
      censor[, i][censor[, i-1] == 1] <- 1
    }
  }


  if (format == "long") {

    statins_long <- reshape(data.frame(statins), varying = 1:total_followup,
                            v.names = "statins", times = seq(start_year, start_year+total_followup-1, by = 1), direction = "long")
    hyper_long <- reshape(data.frame(hyper), varying = 1:total_followup,
                            v.names = "hyper", direction = "long")
    bmi_long <- reshape(data.frame(bmi), varying = 1:total_followup,
                            v.names = "bmi", direction = "long")
    obsdata = data.frame(id = statins_long$id, time = statins_long$time, statins = statins_long$statins, hyper = hyper_long$hyper, bmi = bmi_long$bmi,age = age, sex=sex, y = y)
  }


    if (include_censor && format == "long"){
      obsdata$censor <- as.vector(censor)
    }



  if (include_censor && format == "wide"){
    colnames(censor) <- paste0("censor", seq(start_year, start_year+total_followup-1, by = 1))
    obsdata <- cbind(obsdata, censor)
  }


}

  if(timedep_outcome == TRUE){

    y <- matrix(NA, nrow = n, ncol = total_followup)
    y[, 1] <- rbinom(n, 1, plogis(-2.5 -0.5 * statins[, 1] + 0.25* hyper[, 1] + 0.25 * bmi[, 1]))

    for (i in 2:total_followup) {
    y[,i] <- rbinom(n, 1, plogis(-2.5 -0.5 * statins[, i] + 0.25* hyper[, i] + 0.25 * bmi[, i]))
    y[,i][y[, i-1] == 1] <- 1
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

    obsdata <- data.frame(id = id, age = age, sex = sex,statins = statins, hyper = hyper, bmi = bmi,y = y)

    coln <- paste0(rep(c("statins", "hyper", "bmi","y"), each = total_followup), seq(start_year, start_year+total_followup-1, by = 1))
    colnames(obsdata)[4:ncol(obsdata)] <- coln

    if (format == "long") {
      statins_long <- reshape(data.frame(statins), varying = 1:total_followup,
                              v.names = "statins",times =seq(start_year, start_year+total_followup-1, by = 1),  direction = "long")
      hyper_long <- reshape(data.frame(hyper), varying = 1:total_followup,
                            v.names = "hyper", direction = "long")
      bmi_long <- reshape(data.frame(bmi), varying = 1:total_followup,
                          v.names = "bmi", direction = "long")
      y_long <- reshape(data.frame(y), varying = 1:total_followup,
                          v.names = "y", direction = "long")
      obsdata = data.frame(id = statins_long$id, time = statins_long$time, statins = statins_long$statins,
                          hyper = hyper_long$hyper, bmi = bmi_long$bmi,age = age, sex=sex, y = y_long$y)

    }


      if (include_censor  && format == "long") {
        obsdata$censor <- as.vector(censor)
      }


    if (include_censor && format == "wide") {
      colnames(censor) <- paste0("censor", seq(start_year, start_year+total_followup-1, by = 1))
      obsdata <- cbind(obsdata, censor)
    }


  }
return(obsdata)
}

