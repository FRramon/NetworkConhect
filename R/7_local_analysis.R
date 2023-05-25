################################################################################
# THIS SCRIPT IS DEDICATED TO LOCAL LEVEL ANALYSIS ON BRAIN NETWORK
################################################################################

#' Compute the measure of choice for each node for each patient in a group ("V1", "V2", or "V3)
#'
#' From the dataframe containing all information for one weighting scheme, this function returns a vector
#' of n values by m nodes, where n is the number of participant in a group. And the values are global measure of network topology.
#'
#' @param data dataframe containing all information for one weighting scheme
#' @param v_id chr visit id, like "V1", "V3"
#' @param metric chr the weighting scheme
#' @param eval chr the measure to access the network topology
#' @param threhs_method chr thresh_method is either "threshold" or "density".
#' @param tvalue float either a threshold value or a density value
#' @returns vector of size n x m
#' @examples
#' computeSW(dataFA,"V1","FA","global_eff",0.2)
#' @export
computeSparseMetric <- function(data,
                                v_id,
                                WM_metric = c('FBC','ODI','GFA','Fintra','FA'),
                                eval = c("degree","clust_coeff","strength","betweenness","closeness","eigen","efficiency"),
                                thresh_method,
                                tvalue
                                ){

  WM_metric<-match.arg(WM_metric)
  eval<-match.arg(eval)
  subject_ids <- get_subject_ids(data,v_id)
  RES <- c()
  for(j in 1:length(subject_ids)){
    s_id <- subject_ids[j]
    if(thresh_method == "threshold"){
      g <- makeGraph(data,s_id,v_id,WM_metric,tvalue)
    }
    else if(thresh_method == "density"){
      g0 <- makeGraph(data,s_id,v_id,WM_metric,0)
      g <- sparseThresh(g0,tvalue)
    }
    if (eval =='degree'){
      value <- degree(g)
    } else if (eval == "clust_coeff"){
      # For clustering coeff, all nodes must have at least 2 neighbors.
      # Because there is (n-1) at the denominator
      Isolated = which(degree(g) <= 2)
      g <- delete_vertices(g,Isolated)
      value <- transitivity(g,'weighted')
      names(value) <- V(g)$name
    } else if (eval =='strength'){
      value <- strength(g)
    } else if (eval =='betweenness'){
      value <- betweenness(g,weights = 1/E(g)$weight)
    } else if (eval == 'closeness'){
      value <- closeness(g,weights = 1/E(g)$weight)
    } else if (eval == 'eigen'){
      value <- eigen_centrality(g)$vector
    } else if (eval == 'efficiency'){
      value <- closeness(g,)
    }
    RES<- cbind(RES,value)
  }
  RES
}
#' Kruskal wallis test for non-normal data
K_analysis <- function(Gdata){
  f <- as.formula(paste(colnames(Gdata[1]), "~ Group"))
  kres <- kruskal.test(f ,
                       data = Gdata)
  kres
}

#' Shapiro Wilk test for normality
#'
#' Perform a Shapiro Wilk test to test the normality hypothesis. To further
#' perform parametric test Anova.
#' @param Gdata dataframe two columns : group, y, created by getData
#' @examples
#' isnormal(GdataFA)
#' @export
isnormal <- function(Gdata){
  f <- as.formula(paste(colnames(Gdata[1]), "~ Group"))
  res_aov <- aov(f , data = Gdata)
  st <- shapiro.test(res_aov$residuals)
  if(st[2]>0.05){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}


#' Evaluate group influence every node-level topology measure : Multiple Kruskal-Wallis
#'
#' Perform m Kruskal Wallis on the three groups, where m is the number of nodes. Dependent variable is the vertex-level topology measure.
#' within factor is the group (V1,V2,V3). You can choose a thresholding method : threshold or density, to perform the analysis on thresholded matrix.
#' @param whole_data the data frame containing all data for one weighting scheme
#' @param metric chr the weighting scheme such as "FA", "GFA"
#' @param eval chr the global network topology measure : like clust_coeff, global_eff
#' @param thresh_method chr thresholding procedure, either "threshold" or "density"
#' @param tvalue float threshold value, either a threshold value or a density value
#' @examples
#' eval_on_single_threshold("FA","global_eff",0.2)
#' @export
compute_on_all_nodes <- function(whole_data,metric,eval,thresh_method,tvalue){

  M1 <- computeSparseMetric(whole_data,'V1',metric,eval,thresh_method,tvalue)
  M2 <- computeSparseMetric(whole_data,'V2',metric,eval,thresh_method,tvalue)
  M3 <- computeSparseMetric(whole_data,'V3',metric,eval,thresh_method,tvalue)

  #print(head(M1))
  all_nodes <- c(rownames(M1),rownames(M2),rownames(M3))
  #print(all_nodes)
  counts_nodes <- table(all_nodes)
  #print(counts_nodes)

  all_nodes_once <- unique(all_nodes)
  cols <- c('node','pvalue','mean')
  comparisons <- c('V1 - V2','V1 - V3','V2 - V3')
  comparisons_trunc <-  c('V1 - V2','V1 - V3')
  for(i in 1:length(all_nodes_once)){
    node <- all_nodes_once[i]

    #print(counts_nodes[node][[1]])
    if (counts_nodes[node][[1]] == 3){
      d1 <- as.list(data.frame(M1[node,]))[[1]]
      d2 <- as.list(data.frame(M2[node,]))[[1]]
      d3 <- as.list(data.frame(M3[node,]))[[1]]

      s1 <- get_subject_ids(whole_data,"V1")
      s2 <- get_subject_ids(whole_data,"V2")
      s3 <- get_subject_ids(whole_data,"V3")

      Gdata <- data.frame(
        Y=c(d1,d2,d3),
        Group =factor(rep(c("V1", "V2", "V3"), times=c(length(d1), length(d2), length(d3))))
        )
       # ids = c(s1,s2,s3))



      #Start post hoc  :let i,j in {1,3}, i<j, 1 if Vi>Vj, -1 if Vi > Vj, 0 if not significant
      mV1 <- mean(d1)
      mV2<- mean(d2)
      mV3<- mean(d3)
      D <- dunnTest(Y~Group,Gdata)
      comp <- c(ifelse(mV1>mV2,1,-1),ifelse(mV1>mV3,1,-1),ifelse(mV2>mV3,1,-1))
      is_sign <- D$res$P.adj<0.05
      c1 <- comp*is_sign
      comparisons <- cbind(comparisons,c1)

      comptrunc <- c(ifelse(mV1>mV2,1,-1),ifelse(mV1>mV3,1,-1))
      is_sign_trunc <- D$res$P.adj[c(1,2)]<0.05
      c1_trunc <- comptrunc*is_sign_trunc
      comparisons_trunc <- cbind(comparisons_trunc,c1_trunc)
      #print(is_sign_trunc)
      node_info <- c(node,K_analysis(Gdata)$p.value,mean(mean(d1),mean(d2),mean(d3)))
      cols <- cbind(cols,node_info)

      # LM <- lmerTest::lmer(Y~Group + (1|ids))
      # coeffs <- LM$coefficients
      # pvalues <- unname(coeffs[,5])
      # estimates <- unname(coeffs[,1])
      # # 1 if significant increase, -1 is significant decrease, 0 otherwise
      # comparison_lmer <- (unname(t[,1])*(unname(t[,5]<0.05))!=0)*1
      # comparisons <- c(comparisons,comparison_lmer)
      #
      # node_info <- c(node,K_analysis(Gdata)$p.value,mean(mean(d1),mean(d2),mean(d3)))
      # cols <- cbind(cols,node_info)

      # end post hoc
    }else{print("not in all group")}
  }
  dcols <-t(subset(data.frame(cols),select = -cols))

  dcols <- as.data.frame(dcols)
  rownames(dcols) <- dcols[,1]
  colnames(dcols) <- c('node','pvalue','mean')

  list("M1" = M1,"M2" = M2, "M3" = M3,"Results" = dcols,"Comparison" = comparisons_trunc)#"
}

#' delete je pense
#' @export
boxplot_singlenode <- function(R,node){

  d1 <- as.list(data.frame(R$M1[node,]))[[1]]
  d2 <- as.list(data.frame(R$M2[node,]))[[1]]
  d3 <- as.list(data.frame(R$M3[node,]))[[1]]
  df <- data.frame(name = c(rep('V1',length(d1)),rep('V2',length(d2)),rep('V3',length(d3))),values = c(d1,d2,d3))

  ggbetweenstats(
    data = df,
    x = name,
    y = values,
    type = "nonparametric", #  Kruskal-Wallis
    plot.type = "boxviolin",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    centrality.plotting = FALSE,
    title = paste("Node ",node),
    xlab = "Group",
    bf.message = FALSE
  )
}

#' Reshape data into the main function
#' @export
get_nodal_res <- function(Res,metric,eval){
  test_comp <- Res$R$Comparison
  d <- t(as.data.frame(test_comp))

  colnames(d) <- d[1,]
  d <- d[c(-1),]
  rownames(d) <- Res$R$Results$labelname

  x <- Res$R$Results$labelname
  y <- colnames(d)
  data <- expand.grid(X = x, Y = y)
  data$Z <- as.numeric(d)
  data
}

#' Plot significant result as a heatmap
#'
#' Plot the post hoc result as a heatmap for V1-V2 and V1-V3.
#' light blue represent a significant decrease between V1 and V2/3.
#' regular blue represent non significant result
#' dark blue represent significant increase between V1 and V2/3
#'
#' @param data dataframe containing pvalues and either -1,0,1 to indicate if there is an increase, nothing or a decrease;
#' @param metric chr weighting scheme : "FA", "FBC"...
#' @param eval chr vertex level measure :"clust_coeff", "efficiency"...
#' @param thresh_method chr thresholding procedure, either "threshold" or "density"
#' @param tvalue float threshold value, either a threshold value or a density value
#'
#' @returns save image in Results/metric/Local
#' @export
plotSaveLocalHeatmap <- function(data, metric, eval, thresh_method,tvalue){
  #print(as.character(tvalue))
  path_file <- paste("/home/imabrain/Documents/Brain_networks/vAllWeights/Results/",metric,"/Local/",metric,"_",eval,thresh_method,as.character(tvalue),".png",sep="")
  #print(path_file)
  png(path_file, width = 1920, height = 1080)
  p1 <- ggplot(data, aes(X,Y, fill = Z)) +
    geom_tile(color = "white") +
    theme(axis.text = element_text(angle=90,vjust=0.5,hjust=1))  +
    coord_fixed(ratio = 4) +
    ggtitle(paste(eval,"in",metric,"weighted network at",as.character(tvalue),thresh_method)) +
    xlab("ROIS") +
    ylab("Post-hoc comparison")
  print(p1)
  dev.off()
  print(p1)
}

##################################
# MAIN
##################################


#' Evaluate group influence every node-level topology measure : Multiple Kruskal-Wallis
#'
#' Perform m Kruskal Wallis on the three groups, where m is the number of nodes. Dependent variable is the vertex-level topology measure.
#' within factor is the group (V1,V2,V3). You can choose a thresholding method : threshold or density, to perform the analysis on thresholded matrix.
#' @param whole_data the data frame containing all data for one weighting scheme
#' @param metric chr the weighting scheme such as "FA", "GFA"
#' @param eval chr the global network topology measure : like clust_coeff, global_eff
#' @param thresh_method chr thresholding procedure, either "threshold" or "density"
#' @param tvalue float threshold value, either a threshold value or a density value
#'
#' @returns results, a plot on a heatmap, a plot on the brain atlases
#' @examples
#' main_local_analysis("FA","global_eff",0.2)
#' @export
main_local_analysis <- function(metric,eval,thresh_method,tvalue){

  #LOAD DATA
  print("Step 1: Loading data...")
  data_path <- getDataDir(metric)
  print(data_path)
  whole_data <-read_and_normalize_data(data_path,metric)

  # EXECUTE FUNCTION
  print("Step 2 : Graph creation & computing measure on all node & statistical analysis")
  R <- compute_on_all_nodes(whole_data,metric,eval,thresh_method,tvalue)
  Msign <- subset(R$Results, pvalue <0.05)
  #boxplot_singlenode(R,"53")

  LUT <- getLUT()
  colnames(Msign) <- c("No","pvalue","value")
  Msign <- join(Msign,LUT,by="No")

  colnames(R$Results) <- c("No","pvalue","value")
  R$Results <- join(R$Results,LUT,by="No")
  #print('ok')
  Res <- list("R" = R,"Msign" = Msign)
  dataC <- get_nodal_res(Res,metric,eval)
  # write.csv(dataC,paste(Results/",metric,"/nodal_",eval,thresh_method,as.character(tvalue),".csv",sep=""))
  # write.csv(R$M1,paste(Results/metric,"/V1_nodal_",eval,thresh_method,as.character(tvalue),".csv",sep=""))
  # write.csv(R$M2,paste(Results/metric,"/V2_nodal_",eval,thresh_method,as.character(tvalue),".csv",sep=""))
  # write.csv(R$M3,paste(/Results/metric,"/V3_nodal_",eval,thresh_method,as.character(tvalue),".csv",sep=""))
  plotSaveLocalHeatmap(dataC, metric, eval, thresh_method,tvalue)
  Res
}

