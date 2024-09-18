# library(igraph)
# library(brainGraph)
# library(FSA)
# library(ggstatsplot)
# library(ggplot2)
# library(PMCMRplus)
# library(reshape2)
# library(lme4)
#
# source("/home/imabrain/Documents/Brain_networks/vAllWeights/scripts/1_load_data.R",local=TRUE)
# source("/home/imabrain/Documents/Brain_networks/vAllWeights/scripts/2_graph_construction.R",local=TRUE)
# source("/home/imabrain/Documents/Brain_networks/vAllWeights/scripts/3_compute_metrics.R",local=TRUE)
#
# merge_G_sex_age <- function(){
#   ages <- read_excel("/home/imabrain/Documents/Brain_networks/vAllWeights/data/ages_conhect.xlsx")
#   sex <- read_excel("/home/imabrain/Documents/Brain_networks/vAllWeights/data/Genres_conhect.xlsx")
#
#
#   AAS <- merge(ages,sex, by = "subject_id")
#   AAS
# }
#
# computeSparseMetric <- function(data,
#                                 v_id,
#                                 WM_metric = c('FBC','ODI','GFA','Fintra','FA'),
#                                 eval = c("degree","clust_coeff","strength","betweenness","closeness","eigen","efficiency"),
#                                 wantedDensity,
#                                 hublist
# ){
#
#   WM_metric<-match.arg(WM_metric)
#   eval<-match.arg(eval)
#   subject_ids <- get_subject_ids(data,v_id)
#   RES <- c()
#   for(j in 1:length(subject_ids)){
#     s_id <- subject_ids[j]
#     g0 <- makeGraph(data,s_id,v_id,WM_metric,0)
#     g <- sparseThresh(g0,wantedDensity)
#     if (eval =='degree'){
#       value <- degree(g)
#     } else if (eval == "clust_coeff"){
#       # For clustering coeff, all nodes must have at least 2 neighbors.
#       # Because there is (n-1) at the denominator
#       Isolated = which(degree(g) <= 2)
#       g <- delete_vertices(g,Isolated)
#       value <- transitivity(g,'weighted', v = hublist)
#       names(value) <- V(g)$name
#     } else if (eval =='strength'){
#       value <- strength(g)
#     } else if (eval =='betweenness'){
#       value <- betweenness(g)
#     } else if (eval == 'closeness'){
#       value <- closeness(g)
#     } else if (eval == 'eigen'){
#       value <- eigen_centrality(g)$vector
#     } else if (eval == 'efficiency'){
#       value <- efficiency(g,"nodal")
#     }
#     RES<- cbind(RES,value)
#   }
#   RES
# }
#
# K_analysis <- function(Gdata){
#   f <- as.formula(paste(colnames(Gdata[1]), "~ Group"))
#   kres <- kruskal.test(f ,
#                        data = Gdata)
#   kres
# }
#
# # COMPUTE A GIVEN MEASURE AT A GIVEN DENSITY AND WEIGHTING SCHEME, FOR EVERY NODE
# compute_on_single_node <- function(whole_data,metric,eval,wantedDensity,nodename){
#   M1 <- computeSparseMetric(whole_data,'V1',metric,eval,wantedDensity,nodename)
#   M2 <- computeSparseMetric(whole_data,'V2',metric,eval,wantedDensity,nodename)
#   M3 <- computeSparseMetric(whole_data,'V3',metric,eval,wantedDensity,nodename)
#
#   all_nodes <- c(rownames(M1),rownames(M2),rownames(M3))
#   counts_nodes <- table(all_nodes)
#
#   all_nodes_once <- unique(all_nodes)
#   cols <- c('node','pvalue','mean')
#   comparisons <- c('V1 - V2','V1 - V3','V2 - V3')
#   #for(i in 1:length(all_nodes_once)){
#   index <- which(all_nodes_once == nodename)
#   node <- all_nodes_once[index]
#   print(node)
#   ids1 <- get_subject_ids(whole_data,"V1")
#   ids2 <- get_subject_ids(whole_data,"V2")
#   ids3 <- get_subject_ids(whole_data,"V3")
#
#   #print(counts_nodes[node][[i]])
#   if (counts_nodes[node][[1]] == 3){
#     d1 <- as.list(data.frame(M1[node,]))[[1]]
#     d2 <- as.list(data.frame(M2[node,]))[[1]]
#     d3 <- as.list(data.frame(M3[node,]))[[1]]
#
#     Gdata <- data.frame(
#       Y=c(d1,d2,d3),
#       subject_id = c(ids1,ids2,ids3),
#       Group =factor(rep(c("V1", "V2", "V3"), times=c(length(d1), length(d2), length(d3))))
#     )
#     # Start post hoc  :let i,j in {1,3}, i<j, 1 if Vi>Vj, -1 if Vi > Vj, 0 if not significant
#     mV1 <- mean(d1)
#     mV2<- mean(d2)
#     mV3<- mean(d3)
#     D <- dunnTest(Y~Group,Gdata)
#     comp <- c(ifelse(mV1>mV2,1,-1),ifelse(mV1>mV3,1,-1),ifelse(mV2>mV3,1,-1))
#     is_sign <- D$res$P.adj<0.05
#     c1 <- comp*is_sign
#     comparisons <- cbind(comparisons,c1)
#     node_info <- c(node,K_analysis(Gdata)$p.value,mean(mean(d1),mean(d2),mean(d3)))
#     cols <- cbind(cols,node_info)
#     # end post hoc
#   }else{print("not in all group")}
#   Gdata
# }
# #dcols <-t(subset(data.frame(cols),select = -cols))
#
# #dcols <- as.data.frame(dcols)
# #rownames(dcols) <- dcols[,1]
# #colnames(dcols) <- c('node','pvalue','mean')
#
# #list("M1" = M1,"M2" = M2, "M3" = M3, "Results" = dcols,"Comparison" = comparisons)
# #}
#
#
#
# #ResClust <- main("FBC","betweenness",0.2,hsV1$labels[,1])
# metric <- "GFA"
# data_path <- getDataDir(metric)
# print(data_path)
# whole_data <-read_and_normalize_data(data_path,metric)
#
# # EXECUTE FUNCTION
# R <- compute_on_single_node(whole_data,metric,"efficiency",0.15,"53")
#
#
# AAS <- merge_G_sex_age()
# G2 <- merge(R,AAS, by = "subject_id")
# G2$Group <- as.character(G2$Group)
# #G2$age <- as.numeric(G2$age)
# #G2$Y <- as.numeric(G2$Y)
#
# #subjects <- G2$subject_id
# #groups <- G2$Group
# #y <- G2$Y
# #age <- G2$age
# #sex <- G2$sex
#
# #G3 <- data.frame(
# #  "subject_id" = subjects,
# #  "Y" = y,
# #  "Group" = groups,
# #  "age" = age,
# #  "sex" = sex
# #)
#
# #G4 <- read.table("/home/imabrain/Documents/Brain_networks/vAllWeights/node53_strength.txt")
#
# fit <- lmer(Y ~ Group*sex + (1|subject_id), data = G2)
#
#
#
