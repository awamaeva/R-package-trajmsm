% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_traj.R
\name{predict_traj}
\alias{predict_traj}
\title{Predict trajectory groups for deterministic treatment regimes}
\usage{
predict_traj(
  identifier,
  total_followup,
  treatment,
  time,
  time_values,
  trajmodel
)
}
\arguments{
\item{identifier}{Name of the column of the unique identifier.}

\item{total_followup}{Number of measuring times.}

\item{treatment}{Name of the time-varying treatment.}

\item{time}{Name of the variable time.}

\item{time_values}{Values of the time variable.}

\item{trajmodel}{Trajectory model built with the observed treatment.}
}
\value{
A data.frame with the posterior probabilities.
}
\description{
Function to predict trajectory groups for deterministic treatment regimes
used with gformula and pooled LTMLE.
}
\author{
Awa Diop, Denis Talbot
}
