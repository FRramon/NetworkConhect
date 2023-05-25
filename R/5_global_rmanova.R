


#' Compute the measure of choice for each patient in a group ("V1", "V2", or "V3)
#'
#' From the dataframe containing all information for one weighting scheme, this function returns a vector
#' of n values, where n is the number of participant in a group. And the values are global measure of network topology.
#'
#' @param data dataframe containing all information for one weighting scheme
#' @param v_id chr visit id, like "V1", "V3"
#' @param metric chr the weighting scheme
#' @param eval chr the measure to access the network topology
#' @param threshold float a threshold if the data is computed on a specific threshold. default = 0
#' @examples
#' computeSW(dataFA,"V1","FA","global_eff",0.2)
#' @export
computeSW <- function(data,
                      v_id,
                      WM_metric = c('FBC','ODI','GFA','Fintra','FA'),
                      eval = c('clust_coeff','characteristic_path','global_eff','local_eff','smallworldeness','richcore','strength','betweenness'),
                      thresh_method,
                      tvalue
                      ){

  WM_metric<-match.arg(WM_metric)
  eval<-match.arg(eval)
  subject_ids <- get_subject_ids(data,v_id)
  RES <- vector("numeric", length(subject_ids)) # créer une nouvelle liste vide

  for(j in 1:length(subject_ids)){
    s_id <- subject_ids[j]
    if(thresh_method == "threshold"){
      g <- makeGraph(data,s_id,v_id,WM_metric,tvalue)
    } else if(thresh_method == "density"){
      gT <- makeGraph(data,s_id,v_id,WM_metric,0)
      g <- sparseThresh(gT,tvalue)
      # print("ok")
      # A <- as_adjacency_matrix(g)
      # print("ok")
      # g <- graph_from_adjacency_matrix(A,mode ="undirected",weighted = NULL)
      # print("ok")
    }
    if(eval == 'clust_coeff'){
      value <- transitivity(g,"global")
      #value <- normalized_cluster_coeff(g)
    } else if (eval == 'characteristic_path'){
      #value <- mean_distance(g)
      value <- mean_distance(g,weights = 1/(E(g)$weight))
      #value <- normalized_shortest_path(g)
    } else if (eval == 'global_eff'){
     # value <- global_efficiency(g)
      value <- global_efficiency(g,weights = 1/(E(g)$weight))
    } else if (eval == 'local_eff'){
      #value <- average_local_efficiency(g)
      value <- average_local_efficiency(g,weights = 1/(E(g)$weight))
    } else if (eval == 'smallworldeness'){
      value <- smallworldeness(g)
    } else if (eval == 'richcore'){
      value <- rich_core(g)
    } else if (eval =='strength'){
      value <- mean(strength(g))
    } else if (eval =='betweenness'){
      value <- mean(betweenness(g))
    }
    RES[j] <- value
  }
  RES
}

#' Compute normalized scores
#'
#' Find the best transformation in order that V1 follows a normal distribution,
#' And apply this transformation to V2 and V3.
#'
#' @param Gdata dataframe of patients ids, timepoint (group) and score
#'
#' @returns list of normalized scores for each group
#' @export
NormalizeY <- function(Gdata){
  dataV1 <-subset(Gdata, Group == "V1")
  dataV2 <-subset(Gdata, Group == "V2")
  dataV3 <-subset(Gdata, Group == "V3")

  norm_scoreV1 <- bestNormalize::bestNormalize(dataV1$Y,standardize = TRUE)
  norm_scoreV2 <-predict(norm_scoreV1,dataV2$Y)
  norm_scoreV3 <-predict(norm_scoreV1,dataV3$Y)

  #concatenate
  # Gdata$norm_score <- c(norm_scoreV1$x.t,norm_scoreV2,norm_scoreV3)
  #normvalues <- c(as.list(norm_scoreV1$x.t),as.list(norm_scoreV2),as.list(norm_scoreV3))
  #print(normvalues)

  list("v1" = norm_scoreV1$x.t,"v2" = norm_scoreV2, "v3" = norm_scoreV3,"data" = Gdata)
}

#' Add normalized data column to the dataframe of the results
#'
#' @param Gdata dataframe of patients ids, timepoint (group) and score
#' @param listScores list of normalized scores with respect to V1
#'
#' @returns Gdata with an additional column of normalized data
#' @export
getNormalizeData<-function(Gdata,listScores){
  norm_values <- c(listScores$v1,listScores$v2,listScores$v3)
  Gdata$norm_score <- norm_values
  Gdata
}

#' Group all measures by group, and normalize to normal distribution
#'
#' Compute each measure of network topology for each visit group (V1,V2,V3). Group this into
#' a two columns data frame : group, y. Finally, use bestNormalize to normalize y, to ensure that
#' the normality hypothesis is verified in Anova
#'
#' @param data dataframe containing all information for one weighting scheme
#' @param metric chr the weighting scheme
#' @param eval chr the measure to access the network topology
#' @param threshold float a threshold if the data is computed on a specific threshold. default = 0
#' @examples
#' getData(dataFA,"FA","global_eff",0.2)
#' @export
getData <- function(data,WM_metric,eval,thresh_method,tvalue){
  print("V1")
  gV1 <- computeSW(data,"V1",WM_metric,eval,thresh_method,tvalue)
  gV2 <- computeSW(data,"V2",WM_metric,eval,thresh_method,tvalue)
  gV3 <- computeSW(data,"V3",WM_metric,eval,thresh_method,tvalue)

  idsV1 <- get_subject_ids(data,"V1")
  idsV2 <- get_subject_ids(data,"V2")
  idsV3 <- get_subject_ids(data,"V3")

  Gdata <- data.frame(
    ids = c(idsV1,idsV2,idsV3),
    Y=c(gV1, gV2, gV3),
    Group =factor(rep(c("V1", "V2", "V3"), times=c(length(gV1), length(gV2), length(gV3))))
  )
  #normalize data
 # Gdata$normY <- bestNormalize::bestNormalize(Gdata$Y,standardize = TRUE)$x.t
  #Gdata$normA <- log(Gdata$Y)
 # listScores <- NormalizeY(Gdata)
  #Gdata <- getNormalizeData(Gdata,listScores)
  Gdata
}

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


