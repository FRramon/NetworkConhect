% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_graph_construction.R
\name{makeSparseGraph}
\alias{makeSparseGraph}
\title{Create a igraph object representing a subject at one visit, form one weighting scheme, at a given density}
\usage{
makeSparseGraph(data, s_id, v_id, WM_metric, wantedDensity)
}
\arguments{
\item{data}{dataframe containing all data for one weighting scheme}

\item{s_id}{chr the subject id: like "1-101-MG"...}

\item{v_id}{chr the visit id like "V1"...}

\item{wantedDensity}{the density wanted for the graph}

\item{metric}{chr the weighting scheme : "FA", "FBC", "GFA"...}
}
\value{
a igraph object
}
\description{
Create a weighted undirected graph with igraph. The graph represent a patient (ex : "1-E-101-MG"), at one visit (ex : "V1"),
for one weighting scheme on the edges (ex : "FA"). You need as well to choose a density value. This will truncate the graph by
deleting weakest edges to obtain a graph of the wanted density.
}
