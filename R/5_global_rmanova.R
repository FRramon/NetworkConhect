

#' Shapiro Wilk test for normality
#'
#' Perform a Shapiro Wilk test to test the normality hypothesis. To further
#' perform parametric test Anova.
#' @param Gdata dataframe two columns : group, y, created by getData
#' @examples
#' isnormal(GdataFA)
#' @export
isnormalGlobal <- function(Gdata){
  #Shapiro Wilk test to verify normality
  f <- as.formula("Y ~ Group")
  res_aov <- aov(f , data = Gdata)
  st <- shapiro.test(res_aov$residuals)
  if(st[2]>0.05){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

#' Repeated Measure Anova
#'
#' Perform a RMANOVA on the three groups. Dependent variable is the topology measure.
#' within factor is the group (V1,V2,V3).
#' @param Gdata dataframe two columns : group, y, created by getData
#' @examples
#' RMANOVA_analysis(GdataFA)
#' @export
RMANOVA_analysis <- function(Gdata){
  #repeated measure anova
  res.aov <- anova_test(data = Gdata, dv = norm_score, wid = ids, within = Group)
  #unpaired t comparison with bonferroni correction
  upc <- Gdata %>%
    pairwise_t_test(
      norm_score ~ Group, paired = TRUE,
      p.adjust.method = "bonferroni"
    )
  list("resaov" = res.aov, "upc" = upc)
}


################################################################################
# MAIN
################################################################################


#' Evaluate group influence on one global network topology measure : RMANOVA
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
main_global_analysis <- function(weighting_scheme = c('GFA','FBC','FA','Fintra','ODI'),
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
  # test truncate G for patients that underwent all measurments
  table(G$ids)
  sub_ids <- rownames(data.frame(which(table(G$ids)==3)))
  G <- G[G$ids %in% sub_ids,]
  rownames(G) <- 1:nrow(G)
 # return(G)
  print("Step 3 : Statistical analysis...")
  print(isnormalGlobal(G))
  stats <- RMANOVA_analysis(G)
  print("Step 4 : Results...")
  resaov <- stats$resaov
  print(resaov)
  upc <- stats$upc
  print(upc)

  #par(mfrow = c(2,1))
  #mean_data <- aggregate(score ~ Group,data = G, FUN = mean)
  #print(mean_data)
  p1 <- ggplot(data = G, aes(x = Group,y = norm_score ,group = ids)) +
    theme_light() +
    geom_line(colour = "lightgray") +
    geom_point() +
   # geom_point(data = mean_data,aes(x = Group, y = score,color = Group),size = 3, shape = 19) +
    #geom_line(data = mean_data,aes(x = Group, y = score,color = Group),linetype = "solid") +
    ylab(eval) +
   # stat_summary(fun.y = mean,geom = "line",lwd = 2) +
    ggtitle(paste(eval,"at V1, V2, V3 in",as.character(tvalue),thresh_method))

 upc <- upc %>% add_xy_position(x = "Group")
 #return(upc)
 p2<- ggviolin(G,x = "Group",y = "norm_score",add = "point",fill = 'Group') +
     stat_pvalue_manual(upc) +
     labs(
      subtitle = get_test_label(resaov,detailed = TRUE),
      caption = get_pwc_label(upc)
     ) +
     ylab(eval) +
     ggtitle(paste("RMANOVA on",eval, "at",as.character(tvalue),thresh_method))
 pagePlot <- gridExtra::grid.arrange(p1,p2,ncol =2, nrow = 1)
 print(pagePlot)
 G
}

#main_global_analysis("FA","global_eff","threshold",0.2)


