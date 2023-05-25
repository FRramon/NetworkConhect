# ################################################################################
# # LOAD PACKAGES
# ################################################################################
#
# library(readxl)
# library(igraph)
# library(ggridges)
# library(ggplot2)
# library(shiny)
# library(FSA)
# library(ggstatsplot)
# library(PMCMRplus)
# library(rstatix)
# library(ggpubr)
# library(bestNormalize)
# library(lme4)
# library(reshape2)
# library(dplyr)
# library(ggseg)
# library(ggsegDesterieux)
# library(ggseg3d)
# library(plotly)
# library(plyr)
# library(viridis)
# #library(NetworkConhect)
# library(sjPlot)
# library(glmmTMB)
# options(warn = -1)
# library(lmerTest)
#
# ################################################################################
# #                INDICATE LOCATION OF THE DATA FOLDER                          #
# ################################################################################
# data_path <- "/home/imabrain/Documents/GraphConhect/data"
#
# ################################################################################
# #                           SELECT PARAMETERS                                  #
# ################################################################################
#
# metric <- "FBC"
# thresh_method <- "density"
# threshold <- 0
#
# ################################################################################
# #                            GLOBAL ANALYSIS                                   #
# ################################################################################
#
# ResGlob_g <- main_global_llm(metric, "global_eff", thresh_method,threshold)
# ResGlob_l <- main_global_llm(metric, "characteristic_path", thresh_method,threshold)
# ResGlob_c <- main_global_llm(metric, "clust_coeff", thresh_method,threshold)
# ResGlob_le <- main_global_llm(metric,"local_eff", thresh_method,threshold)
#
# pg<- plot_model(ResGlob_g$fit,show.p = T,show.legend = T,show.data = T,show.values = T,vline.color="red",title = paste("Fixed effects. Y ~ Group + (1|id) ,","global efficiency", "in",metric, "network at",threshold,thresh_method))
# pl<- plot_model(ResGlob_l$fit,show.p = T,show.legend = T,show.data = T,show.values = T,vline.color="red",title = paste("Fixed effects. Y ~ Group + (1|id) ,","characteristic path length", "in",metric, "network at",threshold,thresh_method))
# pc<- plot_model(ResGlob_c$fit,show.p = T,show.legend = T,show.data = T,show.values = T,vline.color="red",title = paste("Fixed effects. Y ~ Group + (1|id) ,","clustering coefficient", "in",metric, "network at",threshold,thresh_method))
# ple<- plot_model(ResGlob_le$fit,show.p = T,show.legend = T,show.data = T,show.values = T,vline.color="red",title = paste("Fixed effects. Y ~ Group + (1|id) ,","local efficiency", "in",metric, "network at",threshold,thresh_method))
#
# multiplot <- gridExtra::grid.arrange(pg,pl,pc,ple,nrow = 2, ncol= 2)
#
#
# ################################################################################
# #                        SAVE RESULTS IN A RESULT FOLDER                       #
# ################################################################################
#
# #dirname(data_path)
#
#
#
