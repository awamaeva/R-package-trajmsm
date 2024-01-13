#' @title Predict trajectory groups for deterministic treatment regimes
#' @description Function to predict trajectory groups for deterministic treatment regimes
#' used with gformula and pooled LTMLE.
#' @name predict_traj
#' @param identifier Name of the column for unique identifier.
#' @param total_follow_up Number of measuring times.
#' @param treatment Name of the time-varying treatment.
#' @param time Name of the measuring times.
#' @param traj_model Trajectory model built with the observed treatment.
#' @import flexmix e1071
#' @export predict_traj
#' @author Awa Diop, Denis Talbot

predict_traj <- function(identifier, total_follow_up, treatment, time, traj_model) {
  data_combn <- bincombinations(total_follow_up)

  # Treatment regime in a long format
  rdata_counter <- reshape(data.frame(data_combn), direction = "long", varying = 1:total_follow_up, sep = "")
  colnames(rdata_counter) <- c(time, treatment, identifier)

  post_probs <- posterior(traj_model, newdata = data.frame(rdata_counter))
  rdata_counter$post_class <- apply(post_probs, 1, which.max)

  data_class <- aggregate(as.formula(paste0("post_class ~ ", identifier)), FUN = unique, data = rdata_counter)

  return(data_class)
}
