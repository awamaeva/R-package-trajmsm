#' @title Split observed data into multiple subsets
#' @description function to split the data into multiple subsets of size s each
#' one corresponding to one time-interval.
#' @name split_data
#' @param obsdata observed data in wide format.
#' @param identifier identifier of individuals.
#' @param ntimes_length number of measuring times per interval.
#' @param total_followup total length of follow-up.
#' @param time name of the time variable.
#' @return \item{all_df}{all subsets, list of time intervals.}
#' @export
#' @author Awa Diop Denis Talbot
#' @examples
#' obsdata = gendata(n = 1000, format = "long", total_followup = 8, seed = 945)
#' years <- 2011:2018
#' res = split_data(obsdata = obsdata, total_followup = 8,
#' ntimes_interval = 6,time = "time", time_values = years,identifier = "id")

split_data <- function(obsdata,total_followup,ntimes_interval,time,time_values, identifier){
total_followup = total_followup
ntimes_interval = ntimes_interval
tot_int <- total_followup -ntimes_interval+1
df_times <- sapply(time_values,function(i)i:(i+ntimes_interval-1))[,1:tot_int]

all_df<- lapply(1:tot_int, function(j){
  all_sub <- obsdata[obsdata[, time] %in% df_times[,j],]
  all_sub$Interv <- j
  all_sub$time2 <- ave(all_sub[,identifier],all_sub[,identifier],FUN = function(x) seq.int(x))
  all_sub$identifier2 <- as.integer(as.character(paste0(all_sub[,identifier], j)))
  return(all_sub)})

return(all_df)
}
