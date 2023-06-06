#' @title Inverse of Probability Weighting
#' @description Compute stabilized and unstabilized with and without censor weights.
#' @name IPW
#' @param id  name of the column for unique identifiant.
#' @param V baseline covariates.
#' @param L time-varying covariates.
#' @param A time-varying treatment.
#' @param K number of measuring times.
#' @param censor name of the censoring variable.
#' @param C logical value TRUE/FALSE to include or not a censoring variable.
#' @param obsdata observed data in wide format.
#' @return \item{IPW}{Stabilized and unstabilized inverse of probabilities with and without censoring}
#' @author Awa Diop, Denis Talbot
#' @export IPW
#' @examples
#' dat = gendatTrajMSM(n=500, Censor=FALSE)
#' V <- c("Age","Sex")
#' L <- c("Hyper","BMI")
#' timevar <- "Time"
#' A <- "Statins"
#' Censor = "C"
#' sw = IPW(numerator = "stabilized",id= "ID", V = c("Age","Sex"),
#'      L = c("Hyper","BMI"),
#'      timevar ="Time",
#'      A=c("Statins"),obsdata = dat)
#' summary(sw)


IPW <- function(numerator = c("stabilized", "unstabilized"),id = id,
                V = V,L = L,A = A, C = FALSE, K = K,
                censor = NULL, time = time, obsdata = obsdata_msm){

if(numerator == "unstabilized" & C == FALSE){
    IPW = IPTW_w(id = id,V = V,L = L,A = A,K = K, time = time, obsdata = obsdata)
}

if(numerator == "stabilized" & C == FALSE){
    IPW = IPTW_sw(id = id,V = V,L = L,A = A,K = K,  time = time, obsdata = obsdata)
}

if(numerator == "unstabilized" & C == TRUE){
    IPW = IPCW_w(id = id,V = V,L = L,A = A,K = K, censor = censor, time = time, obsdata = obsdata)
  }

  if(numerator == "stabilized" & C == TRUE){
    IPW = IPCW_sw(id = id,V = V,L = L,A = A,K  = K, censor = censor, time = time, obsdata = obsdata)
  }
class(IPW) <- "IPW"
return(IPW)
}
