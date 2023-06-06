#' @title ggplotraj
#' @description Use \code{"ggplot2"} to plot trajectory groups produced by the function \code{"buildtraj"} using the observed treatment.
#' @param Rdat an object produced by the function \code{"longtowide"}.
#' @param trt name of the time-varying treatment.
#' @param time name of the variable measurements of time.
#' @param id name of the id variable.
#' @param class name of the trajectory groups.
#' @param FUN specify what statistics to display, by default calculate the mean.
#' @param dpost matrix of the posterior probabilities and the trajectory groups.
#' @param Trajdat merged datasets containing observed data in long format and trajectory groups.
#' @param \dots to add supplementary functions.
#' @return \item{p}{Plot the trajectory groups using the observed treatment.}
#' @author Awa Diop, Denis Talbot
#' @import ggplot2 flexmix
#' @export ggplotraj
#' @examples
#' Obsdata = gendatatraj()
#' Rdat = longtowide(Obsdata = Obsdata, varying = 1:5)
#' head(Rdat)
#' res.traj = buildtraj(Rdat = Rdat, k=3,formula = cbind(A,1-A) ~ time, id="id")
#' dpost = res.traj$dpost
#' head(dpost)
#' ggplotraj(Rdat = Rdat, dpost = dpost, formula = A ~ time + class,
#' trt = "A",time = "time",id="id",class = "class", FUN = "mean")


ggplotraj <- function(Rdat = NULL, dpost = NULL, trajdat = NULL,
                    formula = as.formula(A ~ time + class),
                    trt = "A",time = "time",id="id",class = "class", FUN = "mean",...){

  args <- match.call()

  #Input checks
  if(!is.null(Rdat) & !is.null(trajdat)) stop("Specify only Rdat/dpost or trajdat")
  if(!is.null(dpost) & !is.null(trajdat)) stop("Specify only Rdat/dpost or trajdat")
  if(!("formula" %in% names(args))) stop("Formula not specified")
  if(!("trt" %in% names(args))) stop("Treatment/exposure name not specified")
  if(!("time" %in% names(args))) stop("time name not specified")
  if(!("id" %in% names(args))) stop("id name not specified")
  if(!("class" %in% names(args))) stop("class name not specified")
  if(!("FUN" %in% names(args))) stop("Function not specified")

  if(is.null(trajdat)){

  #Handle variable names with as.symbol()
    trt <- as.symbol(trt)
    time <- as.symbol(time)
    id <- as.symbol(id)
    class <- as.symbol(class)
  #eval(deparse()) to recover variable names
  Rdat$id <- as.integer(as.character(Rdat[,eval(deparse(id))]))
  dpost$id <- as.integer(as.character(dpost[,eval(deparse(id))]))
  Trajdat <- merge(Rdat, dpost, by = eval(deparse(id)))
  #Adherence probabilities
  traj <- aggregate(formula, FUN = FUN, data = Trajdat)

  p <- ggplot(data = traj, aes(x =  !!time, y = !!trt, ...))
  p <- p  + ggtitle("Trajectory groups") + theme_bw()
  p <- p  + geom_point(aes(shape = !!class), size = 3.1) + geom_line(aes(linetype = !!class), size = 1.1)
  p <- p + labs(linetype="Trajectory groups",x="Months post-index",
               y="Adherence probabilities", shape="Trajectory groups")
  }
  if(!is.null(trajdat)){
    #Handle variable names with as.symbol()
    trt <- as.symbol(trt)
    time <- as.symbol(time)
    id <- as.symbol(id)
    class <- as.symbol(class)

    p <- ggplot(data = trajdat, aes(x =  !!time, y = !!trt, ...))
    p <- p  + ggtitle("Trajectory groups") + theme_bw()
    p <- p  + geom_point(aes(shape = !!class), size = 3.1) + geom_line(aes(linetype = !!class), size = 1.1)
    p <- p + labs(linetype="Trajectory groups",x="Months post-index",
                  y="Adherence probabilities", shape="Trajectory groups")
  }
p
}
