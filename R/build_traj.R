#' @title Wrapper for flexmix
#' @description Call the package flexmix to build trajectory groups
#' @name build_traj
#' @param obsdata Data to build trajectory groups in long format.
#' @param number_traj An integer to fix the number of trajectory groups.
#' @param formula Designate the formula to model the longitudinal variable of interest.
#' @param identifier A string to designate the column name for the unique identifier.
#' @param family Designate the type of distribution ("gaussian", "binomial", "poisson", "gamma").
#' @param seed Set a seed for replicability.
#' @param control Object of class FLXcontrol.
#' @param ... Additional arguments passed to the flexmix function.
#' @return A list containing the posterior probability matrix and the fitted trajectory model.
#' @import flexmix
#' @export build_traj
#' @examples
#' obsdata_long = gendata(n = 1000,format = "long", total_followup = 6, seed = 945)
#' formula = as.formula(cbind(statins, 1 - statins) ~ time)
#' restraj = build_traj(obsdata = obsdata_long, number_traj = 3, formula = formula, identifier = "id")

build_traj <- function(obsdata, formula, number_traj, identifier, family = "binomial", seed = 945,
                       control = list(iter.max = 1000, minprior = 0), ...) {
  # Input checks
  if (missing(number_traj)) stop("Number of groups not specified")
  if (missing(identifier)) stop("Identifier not specified")

  # Set seed for replicability
  set.seed(seed)

  # Convert identifier to a symbol if it's a string
  if (is.character(identifier)) {
    identifier <- as.symbol(identifier)
  }

  # Model fitting using flexmix
  res_traj <- flexmix(as.formula(paste0(". ~ .", "|", identifier)), k = number_traj,
                      model = FLXMRglm(formula = formula, family = family),
                      data = obsdata, control = control, ...)

  # Extracting posterior probabilities and trajectory groupings
  data_post <- data.frame(posterior(res_traj)) # Posterior probabilities
  data_post$class <- factor(res_traj@cluster)  # Group trajectories
  data_post[, eval(deparse(identifier))] <- as.integer(as.character(res_traj@group))
  data_post <- unique(data_post)

  return(list(data_post = data_post, traj_model = res_traj))
}
