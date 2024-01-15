#' @title Generate Data Trajectories for MSM
#' @param n Number of observations to generate.
#' @param include_censor Logical, if TRUE, includes censoring.
#' @param format Character, either "long" or "wide" for the format of the output data frame.
#' @export gendata_trajmsm
#' @return A data frame with generated trajectories.
#' @examples
#' gendata_trajmsm(n = 100, include_censor = FALSE, format = "wide")

gendata_trajmsm <- function(n, include_censor = FALSE, format = c("long", "wide"), seed) {
  set.seed(seed)

  # Common variables for all scenarios
  id <- 1:n
  age <- rbinom(n, 1, 0.5)
  sex <- rbinom(n, 1, 0.5)

  # Initialize variables
  statins <- hyper <- bmi <- censor <- matrix(NA, nrow = n, ncol = 3)
  bmi[, 1] <- rbinom(n, 1, plogis(0.5 * age + 0.7 * sex))
  hyper[, 1] <- rbinom(n, 1, plogis(.25 * age + 0.7 * sex + 0.1 * bmi[, 1]))
  statins[, 1] <- rbinom(n, 1, plogis(.75 + 0.4 * age + 0.25* sex - 0.1 * bmi[, 1] - 0.2*hyper[, 1]))
  censor[, 1] <-  rbinom(n, 1, plogis(-2.5 + 0.2 * age + 0.01 * sex + 0.1 * bmi[, 1] - 0.2*hyper[, 1] - 0.5*statins[, 1]))
  # Generate data based on conditions
  for (i in 2:3) {
    bmi[, i] <- rbinom(n, 1, plogis(0.5 * age + 0.7 * sex - 0.65 * statins[, i-1]))
    hyper[, i] <- rbinom(n, 1, plogis(0.5 * age + 0.7 * sex + 0.1 * bmi[, i] - 0.75 * statins[, i-1]))
    statins[, i] <- rbinom(n, 1, plogis(.5 + 0.1 * age + 0.1 * sex - 0.1 * bmi[, i] - 0.2*hyper[, i] + 0.75 * statins[, i-1] ))
    if (include_censor) {
      censor[, i] <- rbinom(n, 1, plogis(-2.5 + 0.02 * age + 0.01 * sex - 0.5 * statins[, i] + 0.1 * hyper[, i] + 0.2 * bmi[, i]))
    }
  }

  # Apply censoring
  if (include_censor) {
    for (i in 2:3) {
      #statins[, i][censor[, i-1] == 1] <- NA
      hyper[, i][censor[, i-1] == 1] <- NA
      bmi[, i][censor[, i-1] == 1] <- NA
      censor[, i][censor[, i-1] == 1] <- 1
    }
  }

  y <- rbinom(n, 1, plogis(0.05 - rowSums(sapply(1:3, function(i) c(-0.25, -0.5, -0.75)[i] * statins[, i] + c(0.2, 0.3, 0.4)[i] * hyper[, i] + c(0.1, 0.2, 0.3)[i] * bmi[, i]), na.rm = TRUE)))

  observed_data <- data.frame(id = id, age = age, sex = sex,
                              statins2011 = statins[, 1], statins2012 = statins[, 2], statins2013 = statins[, 3],
                              hyper2011 = hyper[, 1], hyper2012 = hyper[, 2], hyper2013 = hyper[, 3],
                              bmi2011 = bmi[, 1], bmi2012 = bmi[, 2], bmi2013 = bmi[, 3], y = y)

  if (format == "long") {
    observed_data_long <- reshape(observed_data, varying = list(c("statins2011", "statins2012", "statins2013"),
                                                                c("hyper2011", "hyper2012", "hyper2013"),
                                                                c("bmi2011", "bmi2012", "bmi2013")),
                                  v.names = c("statins", "hyper", "bmi"),
                                  idvar = "id", timevar = "time", times = 2011:2013, direction = "long")
    if (include_censor) {
      observed_data_long$censor <- as.vector(censor)
    }
    return(observed_data_long)
  }

  if (include_censor && format == "wide") {
    observed_data$censor2011 <- censor[, 1]
    observed_data$censor2012 <- censor[, 2]
    observed_data$censor2013 <- censor[, 3]
  }

  return(observed_data)
}

obsdata <- gendata_trajmsm(n = 100, include_censor = TRUE, format = "wide", seed = 345)
