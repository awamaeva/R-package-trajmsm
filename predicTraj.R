#' @title Predict trajectory groups for deterministic treatment regimes
#' @description function to predict trajectory groups for deterministic treatment regimes
#' used with gformula and pooled LTMLE.
#' @name predicTraj.
#' @param id  name of the column for unique identifiant.
#' @param trt name of the time-varying treatment.
#' @param time_name name of the measuring times.
#' @param t number of measuring times.
#' @param name of the id column variable.
#' @param trajmodel trajectory model built with the observed treatment.
#' @import flexmix e1071
#' @author Awa Diop, Denis Talbot


predicTraj <- function(id, t, trt, time_name,trajmodel){
  dat.combn = bincombinations(t)
  #Treatment regime in a long format
  RdatCounter <- reshape(data.frame(dat.combn),direction = "long",varying = 1:t,sep="")
  colnames(RdatCounter) <- c(time_name,trt,id)
  post.probs <- posterior(trajmodel, newdata=data.frame(RdatCounter))
  RdatCounter$postclass<-apply(post.probs,1,which.max)
  dclass <- aggregate(as.formula(paste0("postclass","~", id)),FUN=unique,data= RdatCounter)
  return(dclass)
}
