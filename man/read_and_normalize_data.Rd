% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_load_data.R
\name{read_and_normalize_data}
\alias{read_and_normalize_data}
\title{Create a dataframe with columns "subject_id","visit_id","from","to","weight" for the chosen weighting scheme}
\usage{
read_and_normalize_data(data_path, WM_metric)
}
\arguments{
\item{data_path}{chr indicating the path of the "raw" data}

\item{metric}{a chr for the weighting scheme : "FA", "FBC", "GFA"...}
}
\value{
the newly created dataframe
}
\description{
Create a dataframe with columns "subject_id","visit_id","from","to","weight" for the chosen weighting scheme
}
