library(igraph)
#library(FSA)
#library(ggstatsplot)
#library(ggplot2)
#library(PMCMRplus)
library(reshape2)
library(dplyr)
#library(NetworkConhect)


#' Compute the local measure of choice for each patient in a group ("V1", "V2", or "V3)
#' returns a list of list.(list of each n node local metric for each subject)
#'
#'
#' @param data dataframe containing all information for one weighting scheme
#' @param v_id chr visit id, like "V1", "V3"
#' @param WM_metric chr the weighting scheme
#' @param eval chr the measure to access the network topology
#' @param thresh_method either "threshold" for absolut thresholding, or "density" for density thresholding
#' @param tvalue either the value of the absolute threshold or density threshold
#' @export
computeLocalMetric <- function(data,
                                v_id,
                                WM_metric,
                                eval = c("degree","clust_coeff","strength","betweenness","closeness","eigen","efficiency"),
                                thresh_method,
                                tvalue
                                ){

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
      value <- betweenness(g)
    } else if (eval == 'closeness'){
      value <- closeness(g)
    } else if (eval == 'eigen'){
      value <- eigen_centrality(g)$vector
    } else if (eval == 'efficiency'){
      value <- efficiency(g,"nodal")
    }
    RES<- cbind(RES,value)
  }
  RES
}

#' Compute the local measure of choice for each patient in a group ("V1", "V2", or "V3) and perform the group level statistical analysis
#'
#' Fit a linear mixed model for each region of interest (local metric as dependent variable, group as independant variable, and subject id as random effect to account for repetitive data)
#' Use correction for multiple comparison (bonferronni, fdr)
#'
#' @param df dataframe containing all information for one weighting scheme
#' @param group_mist chr visit id, like "V1", "V3"
#' @param WM_metric chr the weighting scheme
#' @param eval chr the measure to access the network topology
#' @param thresh_method either "threshold" for absolute thresholding, or "density" for density thresholding
#' @param tvalue either the value of the absolute threshold or density threshold
#' @export
main_group_local_metrics_analysis <- function(df,
                                              group_list,
                                              WM_metric,
                                              eval = c("degree","clust_coeff","strength","betweenness","closeness","eigen","efficiency"),
                                              thresh_method,
                                              tvalue){

  ## Compute local metrics for all groups
  L1 <- computeLocalMetric(df,group_list[1],WM_metric,eval,thresh_method,tvalue)
  L2 <- computeLocalMetric(df,group_list[2],WM_metric,eval,thresh_method,tvalue)
  L3 <- computeLocalMetric(df,group_list[3],WM_metric,eval,thresh_method,tvalue)

  # get the labels of the nodes
  all_nodes <- c(rownames(L1),rownames(L2),rownames(L3))
  counts_nodes <- table(all_nodes)
  all_nodes_once <- unique(all_nodes)

  nodes_df <- data.frame(structure = all_nodes_once)


  pvalue_list <- c()
  for(node in all_nodes_once){
    # for each node, create the dataframe structure with correspondant metric, group and subject id
    d1 <- as.list(data.frame(L1[node,]))[[1]]
    d2 <- as.list(data.frame(L2[node,]))[[1]]
    d3 <- as.list(data.frame(L3[node,]))[[1]]

    Gdata <- data.frame(
      subject_id = c(get_subject_ids(df,"V1"),get_subject_ids(df,"V2"),get_subject_ids(df,"V3")),
      Y=c(d1,d2,d3),
      visit_id =factor(rep(c("V1", "V2", "V3"), times=c(length(d1), length(d2), length(d3))))
    )
    # Select relevant columns to avoid duplication
#    df_reduced <- df %>% select(subject_id, visit_id, age, sex, area_from)
    #df_reduced <- df %>% select(subject_id, visit_id)


    # Merge age and sex information into Gdata
    #Gdata <- Gdata %>%
     # left_join(df_reduced, by = c("subject_id", "visit_id"))


    # Gdata <- Gdata %>%
    #   left_join(df,by = c("subject_id","visit_id"))


    # Fit a linear mixed model (random effect as subjecct identificator)
    res.lm <- lmer(Y ~ visit_id + (1|subject_id) ,data = Gdata)
    pvalue_list <- c(pvalue_list,Anova(res.lm)[,"Pr(>Chisq)"])

  }

  # Adjust the n linear models pvalues with FWER method, such as bonferroni (restrictive), of false discovery rate (more permissive)
  padjusted_bonferroni <- p.adjust(pvalue_list,method = 'bonferroni')
  padjusted_fdr <- p.adjust(pvalue_list,method = 'fdr')

  res_dataframe <- data.frame(
    roi = all_nodes_once,
    uncorrected_pvalues = pvalue_list,
    pval_bonferroni = padjusted_bonferroni,
    pval_fdr = padjusted_fdr
  )

  sentence_method = paste("Networks were weighted using ",WM_metric," and thresholded using technique ",thresh_method, " with value ",tvalue,". ",eval," was computed for all nodes. For statistical analysis, one linear mixed models wa used for each node, with node local metric set as the dependent variable, visits as the independent variable, and subject identificators as the random variable. Area of the region was added as a covariate. Finally, multiple comparison correction were done using FDR and/or bonferroni.")

  ## Return only significant regions of interest.
  significant_rois <- res_dataframe[res_dataframe$uncorrected_pvalues < 0.05,]
  list("L1" = L1,"L2" = L2,"L3" = L3,"significant_rois" = significant_rois,"method" = sentence_method)

}

#' This function compute metrics for one node in particular
#' It need the main_group_local_analysis function to have run
#'
#' @param df dataframe containing all information for one weighting scheme
#' @param L1 matrix containing as column as subjects, as rows as nodes. L1 is for V1
#' @param L1 matrix containing as column as subjects, as rows as nodes. L1 is for V1
#' @param L1 matrix containing as column as subjects, as rows as nodes. L1 is for V1
#' @param node str containing the name of the node in the atlas
#' @export
get_one_local_result <- function(df,
                                 L1,
                                 L2,
                                 L3,
                                 node
                                 ){

  d1 <- as.list(data.frame(L1[node,]))[[1]]
  d2 <- as.list(data.frame(L2[node,]))[[1]]
  d3 <- as.list(data.frame(L3[node,]))[[1]]

  Gdata <- data.frame(
    subject_id = c(get_subject_ids(df,"V1"),get_subject_ids(df,"V2"),get_subject_ids(df,"V3")),
    Y=c(d1,d2,d3),
    Group =factor(rep(c("V1", "V2", "V3"), times=c(length(d1), length(d2), length(d3))))
  )

  Gdata

}

# K_analysis <- function(Gdata){
#   f <- as.formula(paste(colnames(Gdata[1]), "~ Group"))
#   kres <- kruskal.test(f ,
#                        data = Gdata)
#   kres
# }
#
# isnormal <- function(Gdata){
#   f <- as.formula(paste(colnames(Gdata[1]), "~ Group"))
#   res_aov <- aov(f , data = Gdata)
#   st <- shapiro.test(res_aov$residuals)
#   if(st[2]>0.05){
#     return(TRUE)
#   }
#   else{
#     return(FALSE)
#   }
# }

# # COMPUTE A GIVEN MEASURE AT A GIVEN DENSITY AND WEIGHTING SCHEME, FOR EVERY NODE
# compute_on_all_nodes <- function(whole_data,metric,eval,thresh_method,tvalue){
#   print(eval)
#   M1 <- computeSparseMetric(whole_data,'V1',metric,eval,thresh_method,tvalue)
#   M2 <- computeSparseMetric(whole_data,'V2',metric,eval,thresh_method,tvalue)
#   M3 <- computeSparseMetric(whole_data,'V3',metric,eval,thresh_method,tvalue)
#
#   print(head(M1))
#   all_nodes <- c(rownames(M1),rownames(M2),rownames(M3))
#   print(all_nodes)
#   counts_nodes <- table(all_nodes)
#   print(counts_nodes)
#
#   all_nodes_once <- unique(all_nodes)
#   cols <- c('node','pvalue','mean')
#   comparisons <- c('V1 - V2','V1 - V3','V2 - V3')
#   comparisons_trunc <-  c('V1 - V2','V1 - V3')
#   for(i in 1:length(all_nodes_once)){
#     node <- all_nodes_once[i]
#
#     print(counts_nodes[node][[1]])
#     if (counts_nodes[node][[1]] == 3){
#       d1 <- as.list(data.frame(M1[node,]))[[1]]
#       d2 <- as.list(data.frame(M2[node,]))[[1]]
#       d3 <- as.list(data.frame(M3[node,]))[[1]]
#
#       Gdata <- data.frame(
#         Y=c(d1,d2,d3),
#         Group =factor(rep(c("V1", "V2", "V3"), times=c(length(d1), length(d2), length(d3))))
#       )
#       # Start post hoc  :let i,j in {1,3}, i<j, 1 if Vi>Vj, -1 if Vi > Vj, 0 if not significant
#       mV1 <- mean(d1)
#       mV2<- mean(d2)
#       mV3<- mean(d3)
#       D <- dunnTest(Y~Group,Gdata)
#       comp <- c(ifelse(mV1>mV2,1,-1),ifelse(mV1>mV3,1,-1),ifelse(mV2>mV3,1,-1))
#       is_sign <- D$res$P.adj<0.05
#       c1 <- comp*is_sign
#       comparisons <- cbind(comparisons,c1)
#
#       comptrunc <- c(ifelse(mV1>mV2,1,-1),ifelse(mV1>mV3,1,-1))
#       is_sign_trunc <- D$res$P.adj[c(1,2)]<0.05
#       c1_trunc <- comptrunc*is_sign_trunc
#       comparisons_trunc <- cbind(comparisons_trunc,c1_trunc)
#       print(is_sign_trunc)
#       node_info <- c(node,K_analysis(Gdata)$p.value,mean(mean(d1),mean(d2),mean(d3)))
#       cols <- cbind(cols,node_info)
#       # end post hoc
#     }else{print("not in all group")}
#   }
#   dcols <-t(subset(data.frame(cols),select = -cols))
#
#   dcols <- as.data.frame(dcols)
#   rownames(dcols) <- dcols[,1]
#   colnames(dcols) <- c('node','pvalue','mean')
#
#   list("M1" = M1,"M2" = M2, "M3" = M3, "Results" = dcols,"Comparison" = comparisons_trunc)
# }
#
#
# boxplot_singlenode <- function(R,node){
#
#   d1 <- as.list(data.frame(R$M1[node,]))[[1]]
#   d2 <- as.list(data.frame(R$M2[node,]))[[1]]
#   d3 <- as.list(data.frame(R$M3[node,]))[[1]]
#   df <- data.frame(name = c(rep('V1',length(d1)),rep('V2',length(d2)),rep('V3',length(d3))),values = c(d1,d2,d3))
#
#   ggbetweenstats(
#     data = df,
#     x = name,
#     y = values,
#     type = "nonparametric", #  Kruskal-Wallis
#     plot.type = "boxviolin",
#     pairwise.comparisons = TRUE,
#     pairwise.display = "significant",
#     centrality.plotting = FALSE,
#     title = paste("Node ",node),
#     xlab = "Group",
#     bf.message = FALSE
#   )
# }
#
# # RESHAPE THE DATA
# get_nodal_res <- function(Res,metric,eval){
#   test_comp <- Res$R$Comparison
#   d <- t(as.data.frame(test_comp))
#
#   colnames(d) <- d[1,]
#   d <- d[c(-1),]
#   rownames(d) <- Res$R$Results$labelname
#
#   x <- Res$R$Results$labelname
#   y <- colnames(d)
#   data <- expand.grid(X = x, Y = y)
#   data$Z <- as.numeric(d)
#   data
# }
#
# # PLOT HEATMAP FUNCTION
# plot_and_save <- function(data, metric, eval, thresh_method,tvalue){
#   print(as.character(tvalue))
#   path_file <- paste("/Volumes/LaCie/derivatives/grouped_results",metric,"/Local/",metric,"_",eval,thresh_method,as.character(tvalue),".png",sep="")
#   print(path_file)
#   png(path_file, width = 1920, height = 1080)
#   p1 <- ggplot(data, aes(X,Y, fill = Z)) +
#     geom_tile(color = "white") +
#     theme(axis.text = element_text(angle=90,vjust=0.5,hjust=1))  +
#     coord_fixed(ratio = 4) +
#     ggtitle(paste(eval,"in",metric,"weighted network at",as.character(tvalue),thresh_method,"threshold")) +
#     xlab("ROIS") +
#     ylab("Post-hoc comparison")
#   print(p1)
#   dev.off()
#   print(p1)
# }

##################################
# MAIN
##################################


# main <- function(metric,eval,thresh_method,tvalue){
#
#   #LOAD DATA
#   # data_path <- getDataDir(metric)
#   # print(data_path)
#   whole_data <-read.csv("/Volumes/LaCie/derivatives/grouped_results/FC_aal.csv")
#   colnames(whole_data) <- c('subject_id','visit_id','from_i','to_j','weight',"from","to")
#   whole_data <- whole_data[c("subject_id","visit_id","from","to","weight")]
#
#
#   whole_data$subject_id <- as.character(whole_data$subject_id)
#   whole_data$visit_id <- as.character(whole_data$visit_id)
#
#   whole_data$visit_id[whole_data$visit_id == "1"] <- "V1"
#   whole_data$visit_id[whole_data$visit_id == "2"] <- "V2"
#   whole_data$visit_id[whole_data$visit_id == "3"] <- "V3"
#
#   # EXECUTE FUNCTION
#   R <- compute_on_all_nodes(whole_data,metric,eval,thresh_method,tvalue)
#   Msign <- subset(R$Results, pvalue <0.05)
#   #boxplot_singlenode(R,"53")
#   #
#   # LUT <- getLUT()
#   # colnames(Msign) <- c("No","pvalue","value")
#   # Msign <- join(Msign,LUT,by="No")
#   #
#   # colnames(R$Results) <- c("No","pvalue","value")
#   # R$Results <- join(R$Results,LUT,by="No")
#   # print('ok')
#   # Res <- list("R" = R,"Msign" = Msign)
#   # dataC <- get_nodal_res(Res,metric,eval)
#   # write.csv(dataC,paste("/home/imabrain/Documents/Brain_networks/vAllWeights/Results/",metric,"/nodal_",eval,thresh_method,as.character(tvalue),".csv",sep=""))
#   # write.csv(R$M1,paste("/home/imabrain/Documents/Brain_networks/vAllWeights/Results/",metric,"/V1_nodal_",eval,thresh_method,as.character(tvalue),".csv",sep=""))
#   # write.csv(R$M2,paste("/home/imabrain/Documents/Brain_networks/vAllWeights/Results/",metric,"/V2_nodal_",eval,thresh_method,as.character(tvalue),".csv",sep=""))
#   # write.csv(R$M3,paste("/home/imabrain/Documents/Brain_networks/vAllWeights/Results/",metric,"/V3_nodal_",eval,thresh_method,as.character(tvalue),".csv",sep=""))
#   # print('ok')
#   # plot_and_save(dataC, metric, eval, thresh_method,tvalue)
#   #
#   # Res
#   R
# }

# R <- main("FBC","degree","threshold",0.1)
