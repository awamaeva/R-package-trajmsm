#' @title Wrapper of flexmix
#' @description Call the package flexmix to construct trajectory groups
#' @name build_traj
#' @param Rdat  sample data to build trajectory groups. These data are in long format.
#'              Each row represent a person time, column A is a binary data.
#' @param J     an integer to choose the number of trajectory groups.
#' @param formula designate the formula to model the longitudinal variable of interest.
#' @param family  designate the type of distribution "gaussian", "binomial", "poisson" and "gamma".
#' @param control object of class FLXcontrol.
#' @param \dots  to add supplementary functions.
#' @return \item{dpost}{Posterior probability.}
#' @return \item{model}{Fitted trajectory model.}
#' @author Awa Diop, Denis Talbot
#' @import flexmix
#' @export build_traj
#' @examples obsdata = gendatatraj()
#' Rdat =longtowide(obsdata = obsdata, varying = 1:5)
#' head(Rdat)
#' res.traj = buildtraj(Rdat = Rdat, k=3,formula = cbind(A,1-A) ~ time, id=id)
#' head(res.traj$dpost)

build_traj <- function(Rdat, formula = as.formula(cbind(A,1-A) ~ time),
              J, id = "id", family = "binomial",
              control = list(iter.max = 1000, minprior = 0), ...){
  #Input checks
  args <- match.call()
  if(!("J" %in% names(args))) stop("Number of groups not specified")
  if(!("id" %in% names(args))) stop("id not specified")

  id <- as.symbol(id)

  res.traj <- flexmix(as.formula(paste0(". ~ .","|",id)), k = J,
                    model = FLXMRglm(formula = formula,
                                     family = family),
                    data = Rdat,control = control,...)
  dpost <- data.frame(posterior(res.traj)) #Posterior probabilities
  dpost$class <- factor(res.traj@cluster)   #Group trajectories
  dpost[,eval(deparse(id))] <- as.integer(as.character(res.traj@group))
  dpost <- unique(dpost)

  return(list(dpost = dpost, model = res.traj))
}
