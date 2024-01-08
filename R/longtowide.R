#' @title Reshape from long to wide format.
#' @description To convert observed data into a wide format for the g-formula and pooled ltmle.
#' @name longtowide
#' @param obsdata the data to reshape in a wide format.
#' @param v.names time varying variables see \code{reshape}.
#' @return \item{long_dat}{Long format data.}
#' @export longtowide
#' @examples
#' widedata = longtowide(obsdata = gendatrajMSM(n=500), idvar = "ID", timevar = "Time")
#' head(widedata)
#' @author Awa Diop, Denis Talbot


longtowide <- function(obsdata,v.names){
  args <- match.call()
  #Input check
  if(!("obsdata" %in% names(args))) stop("Obsdata is missing, with no default")
  if(!("v.names" %in% names(args))) stop("v.names not specified")

  longtowide_dat <- reshape(obsdata,v.names = v.names, direction = "wide")
  return(longtowide_dat)

}
