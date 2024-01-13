#' @title ggplot Trajectory
#' @description Use \code{"ggplot2"} to plot trajectory groups produced by the function \code{"build_traj"} using the observed treatment.
#' @param obsdata Data to plot trajectory groups.
#' @param treatment Name of the time-varying treatment.
#' @param time Name of the variable measurements of time.
#' @param identifier Name of the identifier variable.
#' @param class Name of the trajectory groups.
#' @param fun Specify what statistics to display, by default calculate the mean.
#' @param data_post Matrix of the posterior probabilities and the trajectory groups.
#' @param traj_data Merged datasets containing observed data in long format and trajectory groups.
#' @param \dots Additional arguments to be passed to ggplot functions.
#' @return A ggplot object representing the trajectory groups using the observed treatment.
#' @import ggplot2 flexmix
#' @export ggtraj
#' @examples
#' Obsdata = gendata_trajmsm(n = 1000, include_censor = FALSE, format = "long", seed = 745)
#' Formula = as.formula(cbind(statins, 1 - statins) ~ time)
#' Restraj = build_traj(obsdata = Obsdata, number_traj = 3, formula = Formula, identifier = "id")
#' Datapost = Restraj$data_post
#' head(Datapost)
#' traj_data_long <- merge(Obsdata, Datapost, by = "id")
#'     AggFormula <- as.formula(paste("statins", "~", "time", "+", "class"))
#'     Aggtraj_data <- aggregate(AggFormula, data = traj_data_long, FUN = mean)
#'     Aggtraj_data
#' #Aggtraj_data with labels
#' traj_data_long[ , "traj_group"] <- factor(ifelse(traj_data_long[ , "class"] == "2" ,"Group1" ,
#' ifelse (traj_data_long[ , "class"]== "1" , "Group2" ,"Group3")))
#' AggFormula <- as.formula(paste("statins", "~", "time", "+", "traj_group"))
#' Aggtraj_data <- aggregate(AggFormula, data = traj_data_long, FUN = mean)
#' ggtraj(traj_data =  Aggtraj_data,
#' treatment = "statins",time= "time",identifier="id",class = "traj_group", FUN = mean)

ggtraj <- function(traj_data, treatment, time, identifier, class, FUN = mean, ...) {

  # Input checks
  if(missing(traj_data)) stop("Specify traj_data")
  # Create ggplot
  p <- ggplot(data = traj_data, aes_string(x = time, y = treatment, group = class, color = class, shape = class, linetype = class)) +
    geom_point(size = 3.1) +
    geom_line(size = 1.1) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Observed Treatment Over Time by Trajectory Group",
         x = "Time",
         y = "Adherence Probability",
         color = "Trajectory Group",
         shape = "Trajectory Group",
         linetype = "Trajectory Group") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", hjust = 0.5))
  p <- p + scale_x_continuous(breaks = function(x) unique(as.integer(traj_data[[time]])),
                              labels = function(x) unique(as.integer(traj_data[[time]])))
  return(p)
}

