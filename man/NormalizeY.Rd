% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5_global_linear_model.R
\name{NormalizeY}
\alias{NormalizeY}
\title{Compute normalized scores}
\usage{
NormalizeY(Gdata)
}
\arguments{
\item{Gdata}{dataframe of patients ids, timepoint (group) and score}
}
\value{
list of normalized scores for each group
}
\description{
Find the best transformation in order that V1 follows a normal distribution,
And apply this transformation to V2 and V3.
}
