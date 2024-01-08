#' @title Split observed data into multiple subsets
#' @description function to split the data into multiple subsets of size s each
#' one corresponding to one time-interval.
#' @name split_data
#' @param obsdata observed data in wide format.
#' @param idvar identifiant of individuals.
#' @param s number of measuring times per interval.
#' @param K total length of follow-up.
#' @param timevar name of the time variable.
#' @return \item{all_df}{all subsets, list of time intervals.}
#' @author Awa Diop Denis Talbot
#' @examples res = split_data(obsdata = dat.obs, K = K, s = s, timevar = "time", idvar = "id")

split_data <- function(obsdata,K,s, timevar,idvar){
K = K
s = s
tot_int <- K-s+1
v <- sapply(1:K,function(i)i:(i+s-1))[,1:tot_int]

all_df<- lapply(1:tot_int, function(j){
  all_sub <- obsdata[obsdata[, timevar] %in% v[,j],]
  all_sub$Interv <- j
  all_sub$time2 <- ave(all_sub[,idvar],all_sub[,idvar],FUN = function(x) seq.int(x))
  all_sub$id2 <- as.integer(as.character(paste0(all_sub[,idvar], j)))
  return(all_sub)})

return(all_df)
}
