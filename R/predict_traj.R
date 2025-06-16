#' @title Predict trajectory groups for deterministic treatment regimes
#' @description Function to predict trajectory groups for deterministic treatment regimes
#' used with gformula and pooled LTMLE.
#' @name predict_traj
#' @param identifier Name of the column of the unique identifier.
#' @param total_followup Number of measuring times.
#' @param treatment Name of the time-varying treatment.
#' @param time Name of the variable time.
#' @param time_values Values of the time variable.
#' @param trajmodel Trajectory model built with the observed treatment.
#' @import flexmix e1071
#' @return A data.frame with the posterior probabilities.
#' @export predict_traj
#' @author Awa Diop, Denis Talbot

predict_traj <- function(identifier, total_followup, treatment, time, time_values, trajmodel) {
  data_combn <- bincombinations(total_followup)

  # Treatment regime in a long format
  rdata_counter  <- reshape(data.frame(data_combn), direction = "long", varying = list(1:total_followup),
                       times = time_values, v.names = treatment, timevar = time,
                       idvar = identifier, sep = "")

  post_probs <- posterior(trajmodel, newdata = rdata_counter)

  rdata_counter$post_class <- apply(post_probs, 1, which.max)

  data_class <- aggregate(as.formula(paste0("post_class ~ ", identifier)), FUN = unique, data = rdata_counter)

  return(data_class)
}
