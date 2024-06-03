# Author: Francois Ramon

################################################################################
#                   Main functions for global analysis                        #
################################################################################


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
#' @export
computeSW <- function(data,
                      v_id,
                      WM_metric,
                      eval = c('clust_coeff','characteristic_path','global_eff','local_eff','smallworldeness','richcore','strength','betweenness','edge_weight'),
                      thresh_method,
                      tvalue
){

  eval<-match.arg(eval)
  subject_ids <- get_subject_ids(data,v_id)
  RES <- vector("numeric", length(subject_ids)) # créer une nouvelle liste vide

  for(j in 1:length(subject_ids)){
    s_id <- subject_ids[j]
    if(thresh_method == "threshold"){
      if(WM_metric=="PearsonCorrel"){
        g <- makeGraphFunc(data,s_id,v_id,WM_metric,tvalue)
      }else {
        g <- makeGraph(data,s_id,v_id,WM_metric,tvalue)
      }

    } else if(thresh_method == "density"){
      if(WM_metric=="PearsonCorrel"){
        gT <- makeGraphFunc(data,s_id,v_id,WM_metric,0)
        g <- sparseThresh(gT,tvalue)
      }else {
        gT <- makeGraph(data,s_id,v_id,WM_metric,0)
        g <- sparseThresh(gT,tvalue)

      }

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
    } else if (eval == "edge_weight"){
      value <- mean(E(g)$weight)
    }
    RES[j] <- value
  }
  RES
}

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
#' @export
computeMetric <- function(data,
                          v_id,
                          WM_metric,
                          eval = c('clust_coeff','characteristic_path','global_eff','local_eff','smallworldeness','richcore','strength','betweenness'),
                          thresh_method,
                          tvalue,
                          rsnet = "All"
){
  #LrsNet <- c("Default","SomMot","SalVentAttn","Vis","Limbic","Cont","DorsAttn")
  #rsNet <- LrsNet[5]


  WM_metric<-match.arg(WM_metric)
  eval<-match.arg(eval)
  subject_ids <- get_subject_ids(data,v_id)
  RES <- vector("numeric", length(subject_ids)) # créer une nouvelle liste vide

  for(j in 1:length(subject_ids)){
    s_id <- subject_ids[j]
    if(thresh_method == "threshold"){
      g <- makeGraph(data,s_id,v_id,WM_metric,tvalue)
      if (rsnet != "All"){
        i_DMN <- grep(rsnet,V(g)$name)
        nodes_DMN <- V(g)$name[i_DMN]
        g <- induced_subgraph(
          g,
          i_DMN
        )
      }
    } else if(thresh_method == "density"){
      gT <- makeGraph(data,s_id,v_id,WM_metric,0)
      g <- sparseThresh(gT,tvalue)
      i_DMN <- grep(rsnet,V(g)$name)
      if (rsnet != "All"){
        nodes_DMN <- V(g)$name[i_DMN]
        g <- induced_subgraph(
          g,
          i_DMN)
      }
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
#' @export
getData <- function(data,WM_metric,eval,thresh_method,tvalue){
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
#' Linear Mixed Model
#'
#' Perform a RMANOVA on the three groups. Dependent variable is the topology measure.
#' within factor is the group (V1,V2,V3).
#' @param Gdata dataframe two columns : group, y, created by getData
#' @export
LinearModel_analysis <- function(Gdata){
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

#' plot the data for a good visualization
#' @param G the dataframe
#' @export
plot_dataviz <- function(G,lab){
  p <- ggbetweenstats(
      data = G,
      x = Group,
      y = Y,
      pairwise.comparisons = F,
      pairwise.display = F,
      ylab = lab,
      results.subtitle = F,
      bf.message = F)
  p
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
#' @export
main_global_llm <- function(data_path,weighting_scheme = c('GFA','FBC','FA','Fintra','ODI','PearsonCorrel'),
                                 eval = c('clust_coeff','characteristic_path','global_eff','local_eff','smallworldeness','richcore','strength','betweenness'),
                                 thresh_method,
                                 tvalue
                                 ){

  file_path <- getDataDir(data_path,weighting_scheme)
  print("Step 1 : Loading data...")
  print(file_path)
  whole_data <-read_and_normalize_data(file_path, weighting_scheme)
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
 # plot_dataviz(G)

  print(summary(fit))
  list("fit"=fit,"data"=G)
}



