################################################################################
# LOAD PACKAGES
#1
# library(readxl)
# #2
# library(igraph)
# library(ggridges)
# library(ggplot2)
# #3 & 4
# install.packages("shiny")
# library(shiny)
# #5
# #library(igraph)
# library(brainGraph)
# library(FSA)
# library(ggstatsplot)
# library(PMCMRplus)
# library(rstatix)
# library(ggpubr)
# library(bestNormalize)

#
# devtools::install("NetworkConhect")
# library(NetworkConhect)
# devtools::load_all("/
#library(NetworkConhect)


main <- function(){

  metric <- readline(prompt = "input metric : ")
  global_eval <- readline(prompt = "input global measure : ")
  threshold <- as.numeric(readline(prompt = "input threshold value : "))
  local_eval <- readline(prompt = "input local measure : ")

  ResGlob<- main_global_analysis(metric, global_eval, "threshold",threshold)
  ResLoc<- main_local_analysis(metric,local_eval,"threshold",threshold)


}
#main()

#tdata <- read.table("/home/imabrain/test1.csv")
#fit <- lmerTest::lmer(score ~ Group + (1|ids), data = tdata)
#summary(fit)
