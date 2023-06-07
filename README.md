# NetworkConhect

This R package provides tools for structural brain network analysis\

## Installation

You can install the dev version of NetworkConhect from Github with:

``` r
devtools::install_github("FRramon/NetworkConhect")
```
## Usage

All functions for network analysis are in the package. 
The main function is in the main folder
A data folder containing .xlsx files with corresponding data is needed;

## Analysis Pipeline

### Construction and thresholding

Given the dataset, the first step is to derive a graph from each connectivity matrix. 
* Nodes are derived from brain region automatically parcellated in freesurfer
* Edges are the connection between regions identified in the tractography
* Weight of the edges are either streamline count (the number of fibers between two regions), streamline length (the mean length of the fibers between two regions), or a microstructural metrics, such as fractional anisotropy, neurite density index..\
These steps are shown in figure 1
<p align="center">
  <img width="675" alt="image" src="https://github.com/FRramon/NetworkConhect/assets/109392345/f6264451-3e22-49d2-aaf7-270965ad8d59">
</p>

<font size = "1"> Figure 1 : First steps, graph construction from connectivity data, graph thresholding to ensure comparison in the following steps </font>

### Computing network topology metrics

Each brain network is evaluated on 4 measures, that are shown in figure 2.

<p align="center">
  <img width="649" alt="image" src="https://github.com/FRramon/NetworkConhect/assets/109392345/41c30c73-54b8-4e06-a982-5ed9d7ddde95">
</p>

### Statistical analysis

To observe potential changes in the topology metrics along the treatment, we use a generalized linear mixed effect model. 
$$Y \sim Timepoint + (1|id)$$

