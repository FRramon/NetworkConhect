# Author: Francois Ramon

#' Compute the smallworldenss of a graph
#'
#'Compute the smallworldeness of a graph in the sense of Hagmann 2008.
#' Smallworldeness is defined as the ratio of the normalized clustering by the normalized characteristic path length;
#' Normalized X is defined by dividing X by X' which is the mean of values obtain from 100 random graphs (with same number of edges/nodes and degree/weight distribution)
#'
#' @param g igraph object
#' @export
smallworldeness <- function(g){

  C <- transitivity(g)#,type = "weighted")
  L <- mean_distance(g)
  print(paste("g clust :",as.character(C)))
  print(paste("g dist :",as.character(L)))

  # cr <-list()
  # lr <-list()
  # for(i in 1:10){
  #   deg <- degree(g)
  #   weights <- E(g)$weight
  #   new_weights <- sample(weights)
  #   gR <- sample_degseq(deg,method = 'simple') # Conservation of number of nodes & deg distribution
  #   E(gR)$weight <- new_weights
  #
  #   # Random_elist <- data.frame(as_edgelist(g))
  #   # Random_elist$weight <- new_weights
  #   # gR <- graph_from_data_frame(Random_elist)
  #   cr[i] <- transitivity(gR)
  #   lr[i] <- mean_distance(gR)
  # }
  # df = as.data.frame(unlist(cbind(unlist(cr),unlist(lr))))
  # colnames(df)<-c("clust","path")
  # df$clust[is.nan(df$clust)]<-0
  # df$path[is.nan(df$path)]<-0
  #
  # CR <- mean(df[,1])
  # LR <- mean(df[,2])
  #print(paste("gR clust :",as.character(CR)))
  #print(paste("gR dist :",as.character(LR)))

  #(C/CR)/(L/LR)
  C/L
}



#' Function to perform a random walk on a weighted graph with inverted weights
#' @param g the weighted graph
#' @export
random_walk_length <- function(g, num_walks = 1000, walk_length = 1000) {
  # Inverse the weights of the edges
  E(g)$weight <- 1 / E(g)$weight

  # Get the number of nodes in the graph
  num_nodes <- vcount(g)

  # Generate random starting points
  starting_points <- sample(1:num_nodes, size = num_walks, replace = TRUE)

  # Perform random walks and calculate their lengths
  walk_lengths <- purrr::map_dbl(
    .x = starting_points,
    .f = function(x) {
      walk <- random_walk(g, start = x, steps = walk_length, mode = "all", stuck = "return")
      length(walk) # Get the length of the walk
    }
  )

  # Calculate the mean length of the random walks
  mean_length <- mean(walk_lengths)

  return(mean_length)
}


#' Compute the density of a given graph
#'
#' @param g graph
#' @export
density <- function(g){
  n_max = (vcount(g)*(vcount(g)-1))/2
  n_real = ecount(g)

  n_real/n_max
}


#' Compute the maximal density to use as a threshold
#'
#' It is computed as the minimal density of graphs in the dataset
#'
#' @param g igraph object
#' @param groups list of groups (ex : "V1","V2")
#' @param WM_metric str metric of white matter for edge weighting
#' @export
mindensity <- function(df,groups,WM_metric,rsnet = 'All'){
  L <- c()
  for (group in groups){
    ids = get_subject_ids(df,group)
    for (id in ids){
      if (rsnet != "All"){
        g <- makeGraph(df,id,group,WM_metric,0)
        i_DMN <- grep(rsnet,V(g)$name)
        nodes_DMN <- V(g)$name[i_DMN]
        gi <- induced_subgraph(
          g,
          i_DMN
        )
        L <- c(L,density(gi))
      } else {
      g <- makeGraph(df,id,group,WM_metric,0)
      L <- c(L,density(g))
      }
    }
  }
  L
}


#' Compute the normalized characteristic path lengthof a graph
#'
#'Compute the normalized characteristic path length of a graph.Normalized characteristic path length is defined by dividing characteristic path length
#'by lambda' which is the mean of values obtain from 100 random graphs (with same number of edges/nodes and degree/weight distribution)
#'
#' @param g igraph object
#' @export
normalized_shortest_path <- function(g){
  L <- mean_distance(g)
  lr <-list()
  for(i in 1:10){
    deg <- degree(g)
    gR <- sample_degseq(deg,method = 'simple')
    E(gR)$weight <- sample(E(g)$weight)
    lr[i] <- mean_distance(gR,directed=F)
  }
  df = as.data.frame(unlist(lr))
  colnames(df)<-"path"
  df$path[is.nan(df$path)]<-0
  print(head(df))
  LR <- mean(df[,1])
  L/LR
}

#' Compute the normalized clustering coefficient of a graph
#'
#'Compute the normalized characteristic path length of a graph.Normalized characteristic path length is defined by dividing characteristic path length
#'by lambda' which is the mean of values obtain from 100 random graphs (with same number of edges/nodes and degree/weight distribution)
#'
#' @param g igraph object
#' @export
normalized_cluster_coeff <- function(g){
  C <- transitivity(g)
  cr <-list()
  for(i in 1:10){
    deg <- degree(g)
    gR <- sample_degseq(deg,method = 'simple')
    E(gR)$weight <- sample(E(g)$weight)
    cr[i] <- transitivity(gR)#directed=F)
  }
  df = as.data.frame(unlist(cr))
  colnames(df)<-"path"
  df$path[is.nan(df$path)]<-0
  print(head(df))
  CR <- mean(df[,1])
  C/CR
}


#' Compute the phinorm coefficient for rich club analysis
#'
#'phi norm is defined as the ratio of phi(k) (the density of a graph with nodes of degree >k), divided by the phi_rand(k),
#'which is the same measure for a random graph (with same number of edges/nodes and degree/weight distribution)
#'
#' @param g igraph object
#' @param k int
#' @export
phinorm <- function(g,k){
  set.seed(123)
  phi <- rich_club_coeff(g,k,weighted = TRUE)$phi
  phir <-list()
  for(i in 1:3){
    deg <- degree(g)
    weights <- E(g)$weight
    new_weights <- sample(weights)
    gR <- sample_degseq(deg,method = 'vl') # Conservation of number of nodes & deg distribution
    Random_elist <- data.frame(as_edgelist(g))
    Random_elist$weight <- new_weights
   # print(head(Random_elist))
    gR <- graph_from_data_frame(Random_elist)
    phir[i] <- rich_club_coeff(gR,k,weighted = TRUE)$phi
   # print(phir[i])

  }
  df = as.data.frame(unlist(phir))
  colnames(df)<- "phi"
  df$phi[is.nan(df$phi)]<-0

  phiR <- mean(df[,1])

  phi/phiR
}

#' Identify modules of a graph and membership of each node
#'
#'Module detection is done with Louvain algorithm, for cluster detection.
#'
#' @param g igraph object
#' @export
module_detection <- function(g){
	gamma <- seq(0.25,2,0.025)
	nc <- vector("numeric",length(gamma))
	for(i in 1:length(gamma)){

		gLouvain <- cluster_louvain(g,resolution = gamma[i])
		print(length(gLouvain))
		nc[i] <- length(gLouvain)
	}
		#nCommunities <- sizes(gLouvain)
	#plot(gLouvain,g)
	#print(nCommunities)
	plot(gamma,nc)
}

#' to delete
hub_detectionVisit <- function(g1,g2,g3){
	StrengthG1 <- strength(g1)
	threshS1 <- mean(StrengthG1) + sd(StrengthG1)
	StrengthG2 <- strength(g2)
	threshS2 <- mean(StrengthG2) + sd(StrengthG2)
	StrengthG3 <- strength(g3)
	threshS3 <- mean(StrengthG3) + sd(StrengthG3)
	colors1 <- ifelse(sort(StrengthG1,decreasing = T)>threshS1, "yellow", "lightyellow")
	colors2 <- ifelse(sort(StrengthG2,decreasing = T)>threshS2, "blue", "lightblue")
	colors3 <- ifelse(sort(StrengthG3,decreasing = T)>threshS3, "red", "coral")

	par(mfrow = c(3,1))
	barplot(sort(StrengthG1,decreasing = T),col = colors1,horiz = F, las = 2)
	abline(h=threshS1)
	barplot(sort(StrengthG2,decreasing = T),col = colors2,horiz =F,las = 2)
  barplot(sort(StrengthG3,decreasing = T),col = colors3,horiz = F,las = 2)
  list('hubsV10' = which(StrengthG1>threshS1),'hubsV2' = which(StrengthG2>threshS2),'hubsV3' = which(StrengthG3>threshS3))
}

#' Hub detection on a graph using different methods
#'
#' Nodes are defined as Hubs if they verify the inequality :
#'             X > mean(X) + sd(X)
#' Where X can be node strength, betweenness centrality, or closeness centrality
#' Provide a horizontal histogram plot
#'
#' @param g igraph object
#' @param method chr method for hub detection "closeness","strength","betweenness"
#' @export
hub_detectionh <- function(g,method,data_path){
  if(method == 'strength'){
    distrib <- strength(g)
  } else if (method == 'betweenness'){
    distrib <-betweenness(g,weights = 1/E(g)$weight)
  } else if (method == 'closeness'){
    distrib <- closeness(g,weights = 1/E(g)$weight)
  } else if (method == 'degree'){
    distrib <- degree(g)
  }
  threshold <- mean(distrib) + sd(distrib)
  colors <- ifelse(sort(distrib,decreasing = F)[120:166]>threshold, "yellow", "lightyellow")

  labels <- as.data.frame(rownames(as.data.frame(sort(distrib,decreasing = F))))
  colnames(labels) <- "No"
  LUT <- getLUT(data_path)
  lbls <- merge(labels,LUT,by="No")
  index <- match(labels$No,lbls$No)
  lbls <- lbls[index,]

  p <- barplot(sort(distrib,decreasing = F)[121:166],axisnames=TRUE,col = colors,horiz = T, las = 2,xlab = method,xlim = c(0,max(distrib)),names = lbls[,2][121:166])
  abline(v=threshold)

  #HUBS def
  hlabels <- as.data.frame(rownames(as.data.frame(which(distrib>threshold))))
  colnames(hlabels) <- "No"
  hlbls <- merge(hlabels,LUT,by="No")

  list("plot" = p,"labels" = hlbls)
}


#' Hub detection on a graph using different methods
#'
#' Nodes are defined as Hubs if they verify the inequality :
#'             X > mean(X) + sd(X)
#' Where X can be node strength, betweenness centrality, or closeness centrality
#' Provide a vertical histogram plot
#'
#' @param g igraph object
#' @param method chr method for hub detection "closeness","strength","betweenness"
#' @export
hub_detectionv <- function(g,method){
  if(method == 'strength'){
    distrib <- strength(g)
  } else if (method == 'betweenness'){
    distrib <-betweenness(g)
  } else if (method == 'closeness'){
    distrib <- closeness(g)
  } else if (method == 'degree'){
    distrib <- degree(g)
  }
  threshold <- mean(distrib) + sd(distrib)
  colors <- ifelse(sort(distrib,decreasing = T)>threshold, "yellow", "lightyellow")

  labels <- as.data.frame(rownames(as.data.frame(sort(distrib,decreasing = T))))
  colnames(labels) <- "No"
  LUT <- getLUT()
  lbls <- merge(labels,LUT,by="No")
  index <- match(labels$No,lbls$No)
  lbls <- lbls[index,]

  #png("/home/imabrain/Documents/test_s.png",height = 1080,width = 1920)

  p <- barplot(sort(distrib,decreasing = T)[1:45],axisnames=TRUE,col = colors,horiz = F, las = 2,ylab = method,names = lbls[,2][1:45])
  abline(h=threshold)

  #HUBS def
  hlabels <- as.data.frame(rownames(as.data.frame(which(distrib>threshold))))
  colnames(hlabels) <- "No"
  hlbls <- merge(hlabels,LUT,by="No")

  list("plot" = p,"labels" = hlbls)
}


#' Hub detection on a graph : three methods
#'
#' Nodes are defined as Hubs if they verify the inequality :
#'             X > mean(X) + sd(X)
#' Where X can be node strength, betweenness centrality, or closeness centrality.
#' The difference with the other hub detection function is that it is less restrictive. We consider a node as a hub
#' if it fullfill the inequality for at least one method. Here the function has been created to work on group average data (three graphs in input,
#' corresponding to each average graph for a visit)
#'
#' @param metric chr weighting scheme like "FA", "GFA"
#' @param R is a list of graph containing G1, G2, G3 (one graph for each visit).
#' @export
combine_hubs <- function(metric,R){
  metric <- "GFA"
  data_path <- getDataDir(metric)

  par(mfrow = c(3,1))#,mar = default_margin)
  hsV1 <- hub_detectionv(R$G1,"strength")
  hbV1 <- hub_detectionv(R$G1,"betweenness")
  hcV1 <- hub_detectionv(R$G1,"closeness")
  mtext(paste("Hubs in V1 -",metric),side = 3, line = -3, outer = TRUE)

  hsV2 <- hub_detectionv(R$G2,"strength")
  hbV2 <- hub_detectionv(R$G2,"betweenness")
  hcV2 <- hub_detectionv(R$G2,"closeness")
  mtext(paste("Hubs in V2 -",metric),side = 3, line =-3, outer = TRUE)

  hsV3 <- hub_detectionv(R$G3,"strength")
  hbV3 <- hub_detectionv(R$G3,"betweenness")
  hcV3 <- hub_detectionv(R$G3,"closeness")
  mtext(paste("Hubs in V3 -",metric),side = 3, line = -3, outer = TRUE)

  HubsV1 <- unique(c(hsV1$labels[,2],hbV1$labels[,2],hcV1$labels[,2]))
  HubsV2 <- unique(c(hsV2$labels[,2],hbV2$labels[,2],hcV2$labels[,2]))
  HubsV3 <- unique(c(hsV3$labels[,2],hbV3$labels[,2],hcV3$labels[,2]))

  l1 <- data.frame(h1)
  l2 <- data.frame(h2)
  l3 <- data.frame(h3)

  print(rownames(l1))
  print(rownames(l2))
  print(rownames(l3))
  u <- unique(c(rownames(l1),rownames(l2),rownames(l3)))


  list("hubsV1" = HubsV1,"hubsV2" = HubsV2,"hubsV3" = HubsV3)
}



