#Author Francois Ramon

#' Idea : Visualisation tool on a Desterieux atlas using ggseg or ggseg3d.
#' Would just need a table of values by zones (like degree of each node, col1 = region_label, col2 = value)
#'
#'

library(ggseg)
library(ggseg3d)
library(dplyr)
library(ggsegDesterieux)
library(ggplot2)

cortical_labels <- brain_labels(desterieux)
subcortical_labels <- brain_labels(aseg)

#Plot Desterieux (freesurfer) cortical representation in 2d with labels
#' import
plot_cortdesterieux2d <- function(){
  plot(desterieux) +
    theme(legend.position = "bottom",
         legend.text = element_text(size = 7)) +
    guides(fill = guide_legend(ncol = 6))

}

#' import
plot_cortdesterieux3d <- function(){
  ggseg3d(atlas = desterieux_3d) %>%
    pan_camera("right lateral")

}




dataOnBrain = data.frame(
  label = c("rh_S_front_inf", "lh_S_front_inf"),
  p = sample(seq(0,.5,.001), 2),
  stringsAsFactors = FALSE)

#Make the plot of the cortical region colored by there value of p (which can be either degree, betweenness etc etc)

dataOnBrain %>%
  brain_join(desterieux, by = "label") %>%
  ggplot() +
  geom_sf(aes(fill = p))
