library(lme4)
#' Linear Mixed Model
#'
#' Perform a RMANOVA on the three groups. Dependent variable is the topology measure.
#' within factor is the group (V1,V2,V3).
#' @param Gdata dataframe two columns : group, y, created by getData
#' @examples
#' RMANOVA_analysis(GdataFA)
#' @export
LinearModel_analysis <- function(Gdata){
  #repeated measure anova
#  print(colnames(Gdata))
  llm  <- lme4::glmer(Y ~ Group + (1|ids),data = Gdata,family = Gamma(link = "log"))
  llm
}

plot_results_lmer <- function(fit){

  ranef_fit <- ranef(fit)$ids[, 1]
  #assumptions
  par(mfrow = c(3,3))
  # Check linearity
  plot(fitted(fit), resid(fit),main ="Linearity plot")
  abline(h = 0, lty = 2)

  # Check independence
  acf(resid(fit))

  # Check homoscedasticity
  plot(fitted(fit), abs(resid(fit)),main = "Homoscedasticity plot")

  #Check normmliy of the residuals
  hist(resid(fit))
  qqnorm(resid(fit))
  qqline(resid(fit))
  print(shapiro.test(resid(fit)))

  # Check the normality of the random effects
  hist(ranef_fit)
  qqnorm(ranef_fit)
  qqline(ranef_fit)
  print(shapiro.test(ranef_fit))


  p <- plot_model(fit)#, type = "est",show.p = T, show.value = T)
  t <- tab_model(fit)
  list("plot" = p)#
 # p
}



check_assumptions <- function(Gtest){
  par(mfrow = c(3,2),mar = c(1,1,1,1))
  hist(Gtest$Y, main = "Distribution, all groups",sub = as.character(shapiro.test(Gtest$Y)$p))
  hist(subset(Gtest,Group =="V1")$Y, main = "Distribution in group V1",sub = as.character(shapiro.test(subset(Gtest,Group =="V1")$Y)$p))
  hist(subset(Gtest,Group =="V2")$Y, main = "Distribution in group V2",sub = as.character(shapiro.test(subset(Gtest,Group =="V2")$Y)$p))
  hist(subset(Gtest,Group =="V3")$Y, main = "Distribution in group V3",sub = as.character(shapiro.test(subset(Gtest,Group =="V3")$Y)$p))

  fit <- glmer(Y ~ Group + (1|ids),data = Gtest,family = Gamma(link = "log"))

  hist(resid(fit), main = "distribution of the residuals",sub = as.character(shapiro.test(resid(fit))$p))
  hist(ranef(fit)$ids[, 1],main = "distribution of the random factor", sub = as.character(shapiro.test(ranef(fit)$ids[, 1])$p))
}

################################################################################
# MAIN
################################################################################

#' Evaluate group influence on one global network topology measure : Linear model
#'
#' Perform a RMANOVA on the three groups. Dependent variable is the topology measure.
#' within factor is the group (V1,V2,V3). You can choose a threshold to perform the analysis on thresholded matrix. default is 0
#' @param weighting_scheme chr the weighting scheme such as "FA", "GFA"
#' @param eval chr the global network topology measure : like clust_coeff, global_eff
#' @param thresh_method chr the method of threshoding. Whether "threshold" or "density"
#' @param threshold float the threshold/density to be applied on each network in the analysis.
#' @examples
#' eval_on_single_threshold("FA","global_eff",0.2)
#' @export
main_global_llm <- function(weighting_scheme = c('GFA','FBC','FA','Fintra','ODI'),
                                 eval = c('clust_coeff','characteristic_path','global_eff','local_eff','smallworldeness','richcore','strength','betweenness'),
                                 thresh_method,
                                 tvalue
                                 ){

  data_path <- getDataDir(weighting_scheme)
  print("Step 1 : Loading data...")
  print(data_path)
  whole_data <-read_and_normalize_data(data_path, weighting_scheme)
  print("Step 2 : Graph creation & computing measures...")
  G <- getData(whole_data, weighting_scheme,eval,thresh_method,tvalue)
  check_assumptions(G)
  print(paste("./",eval,weighting_scheme,".csv"))
  write.csv(G,paste("./",eval,weighting_scheme,".csv",sep = ""))
  #return(G)
  print("Step 3 : Statistical analysis...")
  #print(isnormalGlobal(G))
  fit <- LinearModel_analysis(G)
  print("Step 4 : Results...")
  plot_model(fit,show.intercept = T,show.p = T,show.legend = T,show.data = T,show.values = T,vline.color="red")

  plot_results_lmer(fit)

  print(summary(fit))
  list("fit"=fit,"data"=G)
}



