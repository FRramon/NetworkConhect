#Author : Francois Ramon

###############################################################################################
# CREATE GRAPH, EDGELIST,ADJACENCY MATRIX FROM 1_load_data.R
###############################################################################################

#' Create a igraph object representing a subject at one visit, form one weighting scheme, at a given threshold.
#'
#' Create a weighted undirected graph with igraph. The graph represent a patient (ex : "1-E-101-MG"), at one visit (ex : "V1"),
#' for one weighting scheme on the edges (ex : "FA"). You can as well choose a threshold value. This will truncate the graph by
#' deleting the edges that have a value below the threshold. Default is 0
#'
#' @param data dataframe containing all data for one weighting scheme
#' @param s_id chr the subject id: like "1-101-MG"...
#' @param v_id chr the visit id like "V1"...
#' @param metric  chr the weighting scheme : "FA", "FBC", "GFA"...
#' @param threshold float a threshold value, default = 0
#' @returns a igraph object
#' @export
makeGraph <- function(data,s_id,v_id,WM_metric=c('FBC','ODI','GFA','FA','Fintra'),threshold = 0){
	WM_metric<-match.arg(WM_metric)
	df <- subset(data, subject_id == s_id & visit_id == v_id)
  edgelist <- df[, c('from','to','weight')]
	# edgelist<-data.frame(df[,3],df[,4],df[,5])
	# colnames(edgelist) <- c('from','to','weight')
	elist<-subset(edgelist,weight>threshold)
	#finv <- function(x) 1/x
	#elist$weight <- as.numeric(lapply(elist$weight,finv))
	g<-graph_from_data_frame(elist,directed=FALSE)
  #print(is.weighted(g))
	g
}

#' Create a edge list  representing a subject at one visit, form one weighting scheme, at a given threshold.
#'
#' Create a data frame of three columns : from, to, weight. The edge list represent a patient (ex : "1-E-101-MG"), at one visit (ex : "V1"),
#' for one weighting scheme on the edges (ex : "FA"). You can as well choose a threshold value. This will truncate the graph by
#' deleting the edges that have a value below the threshold. Default is 0
#'
#' @param data dataframe containing all data for one weighting scheme
#' @param s_id chr the subject id: like "1-101-MG"...
#' @param v_id chr the visit id like "V1"...
#' @param metric  chr the weighting scheme : "FA", "FBC", "GFA"...
#' @param threshold float a threshold value, default = 0
#' @returns dataframe
#' @export
makeEdgelist <- function(data,s_id,v_id,WM_metric=c('FBC','ODI','GFA','FA','Fintra'),threshold=0){
	WM_metric<-match.arg(WM_metric)
	df <- subset(data, subject_id == s_id & visit_id == v_id)
	edgelist<-data.frame(df[,3],df[,4],df[,5])
	colnames(edgelist) <- c('from','to','weight')
	elist<-subset(edgelist,weight>threshold)
	elist
}

#' Create a igraph object at a specific wanted density
#'
#' Graph density is defined as the number of edges divided by the possible number of edges (defined as 0.5*nn(nn-1) where nn is the number of nodes).
#' The function calculates the number of edges that is necessary to delete to obtain the wantedDensity, delete the lowest weighted number of edges, and returns
#' a newly created graph, that has the wanted density.
#'
#' @param g igraph object
#' @param wantedDensity the density wanted for the new graph

#' @returns a igraph object at density 'wantedDensity'
#' @export
sparseThresh <- function(g, wantedDensity){
  adj_original <- as_adjacency_matrix(g,attr = "weight",sparse = FALSE )
  adj_zero <- adj_original
  adj_zero[] <- 0
  #print(adj_zero)
  df<-get.data.frame(g)
  n <- vcount(g)
  en <- round((n^2 - n) * wantedDensity/ 2) # Calcul du nombre d'arêtes à conserver pour obtenir la
  # "sparsity" voulue.
  spars <- 2*ecount(g)/(n*n-n) # Calcul de la "sparsity" originale
  df_ordered <- df[order(df$weight,decreasing=TRUE),] # classement des valeurs des arêtes
  #dans l'ordre décroissant.
  to_keep <- df_ordered[1:en,]

  WS <- adj_zero

  for (i in 1:nrow(to_keep)) {
    from <- to_keep$from[i]
    to <- to_keep$to[i]
    weight <- to_keep$weight[i]
    WS[from, to] <- weight
    WS[to, from] <- weight # if the graph is undirected
  }
  gS1 <- graph.adjacency(WS, mode="undirected", weighted=TRUE)

  isolated_nodes <- which(degree(gS1) == 0)
  gS <- delete.vertices(gS1,isolated_nodes)
  #print(vcount(gS1))
  #print(vcount(gS))

  nS <- vcount(gS)
  Nspars <- 2*ecount(gS)/(nS*nS-nS)
  #print(paste("New sparsity ",Nspars))
  gS
}

#' Create a igraph object representing a subject at one visit, form one weighting scheme, at a given density
#'
#' Create a weighted undirected graph with igraph. The graph represent a patient (ex : "1-E-101-MG"), at one visit (ex : "V1"),
#' for one weighting scheme on the edges (ex : "FA"). You need as well to choose a density value. This will truncate the graph by
#' deleting weakest edges to obtain a graph of the wanted density.
#'
#' @param data dataframe containing all data for one weighting scheme
#' @param s_id chr the subject id: like "1-101-MG"...
#' @param v_id chr the visit id like "V1"...
#' @param metric  chr the weighting scheme : "FA", "FBC", "GFA"...
#' @param wantedDensity the density wanted for the graph
#' @returns a igraph object
#' @export
makeSparseGraph <- function(data,s_id,v_id,WM_metric=c('FBC','ODI','GFA','FA','Fintra'),wantedDensity){
  WM_metric<-match.arg(WM_metric)
  df <- subset(data, subject_id == s_id & visit_id == v_id)
  edgelist<-data.frame(df[,3],df[,4],df[,5])
  colnames(edgelist) <- c('from','to','weight')
  g<-graph_from_data_frame(edgelist,directed=FALSE)
  gS <- sparseThresh(g,wantedDensity)
  gS
}

#' Return the minimum graph density across all participants and all visits
#'
#' Compute all graph densities and returns the mininum
#'
#' @param data dataframe containing all data for one weighting scheme
#' @param metric  chr the weighting scheme : "FA", "FBC", "GFA"...
#' @returns float group minimum density
#' @export
getMeanSparsity <- function(data,WM_metric = c('FBC','ODI','GFA','FA','Fintra')){
  timePoint <- c('V1','V2','V3')
  Lspars <- c()
  for(v_id in timePoint){
    subjectIds <- get_subject_ids(data,v_id)
    for(s_id in subjectIds){
      g <- makeGraph(data,s_id,v_id,WM_metric,0)
      spars <- 2*ecount(g)/(vcount(g)*vcount(g)-vcount(g))
      Lspars <- c(Lspars,spars)
    }
  }
  min(Lspars)
}

#' Plot the adjacency matric of a graph
#'
#' This function plot the adjacency matrix of a graph as a heatmap.
#'
#' @param g igraph object
#' @param metric  chr the weighting scheme : "FA", "FBC", "GFA"...
#' @export
plot_adjacency <- function(g,WM_metric=c('FBC','ODI','GFA','FA','Fintra')){
	WM_metric<-match.arg(WM_metric)
	adj <- as_adjacency_matrix(g,attr="weight",sparse = FALSE)
	n<-vcount(g)
	palf <- colorRampPalette(c('blue','yellow','red'))
	if(WM_metric=='FBC'){
	  heatmap(log(adj[,n:1]), Rowv = NA, Colv = NA,scale="none",col = palf(200))
	}else{
	  heatmap(adj[,n:1],Rowv = NA,Colv=NA,scale="none",col = palf(200))
	}
}

#' Plot the edges weight distribution for one participant : FBC
#'
#'TO DELETE OR UPDATE
#' Plot the fiber count edges weight distribution on a log scale for one participant.
#'
#' @param data edge list
#' @export
plot_weight_distribution <- function(data){
	mybreaks <- seq(0,max(data$weight,1000))
	h<-hist(data$weight,breaks = mybreaks,xlab = "Nombre de fibres",ylab = "Nombre de faisceaux",main = "Faisceaux cerveau entier",plot= FALSE)
	plot(h$count,main = "Distribution of fiber counts : Log scaled", xlab = "Log(fiber number)",ylab = "count",log = "x")
}

#' Plot the edges weight distribution for one participant : FBC
#'
#'TO DELETE OR UPDATE
#' Plot the fiber count edges weight distribution on a log scale for one participant. from a igraph object
#'
#' @param g igraph object
#' @param nBreaks number of breaks for the histogram
#' @export
plot_weight_distribution_graph <- function(g,nBreaks){
  w <- data.frame(E(g)$weight)[,1]
  mybreaks <- seq(min(w),max(w)+1000,(max(w)+ 1000 - min(w))/nBreaks)
  mybreaksplot <- seq(min(w),max(w)+1000,(max(w)+ 1000 - min(w))/(nBreaks-1))
  h<-hist(w,breaks = mybreaks,xlab = "Nombre de fibres",ylab = "Nombre de faisceaux",main = "Faisceaux cerveau entier",plot= FALSE)
  plot(mybreaksplot,h$count,main = "Distribution of fiber counts : Log scaled", xlab = "fiber number",ylab = "count",type = 'h',lwd = 10,lend =2,log = "y")
}

#' Plot the edges weight distribution for one participant : FA
#'
#'TO DELETE OR UPDATE
#' Plot the fractional anistropy edges weight distribution for one participant.from a graph
#'
#' @param g igraph object
#' @param nBreaks number of breaks for the histogram
#' @export
plot_weight_distribution_graph_fa <- function(g,v_id,nBreaks){
  w <- data.frame(E(g)$weight)[,1]
  mybreaks <- seq(min(w),max(w),(max(w) - min(w))/nBreaks)
  mybreaksplot <- seq(min(w),max(w),(max(w) - min(w))/(nBreaks-1))
  h<-hist(w,breaks = mybreaks,xlab = "Nombre de fibres",ylab = "Nombre de faisceaux",main = "Faisceaux cerveau entier",plot= FALSE)
  plot(mybreaksplot,h$count,main = paste("Distribution of edge weights - FA -",v_id), xlab = "FA value",ylab = "Count",type = 'h',lwd = 10,lend =2)
}

#' Plot the edges weight distribution for each participant of a group. FA weighted
#'
#'
#' Plot the fractional anistropy edges weight distribution each participants of a group
#' (V1, V2, V3). Plot it as a ridge plot. Note : You can plot the three ridges side to side
#' by calling this function three times and using gridExtra::grid.arrange.
#'
#' @param dataFA dataframe containing all the FA weighted data
#' @param v_id chr visit id like "V1", "V2"...
#' @export
plot_ridge_distribution <- function(dataFA,v_id){
  idsV1 <- get_subject_ids(dataFA,v_id)
  graph_V1 <- lapply( idsV1, FUN =  function(x) makeGraph(dataFA,x,v_id,"FA",0))
  edges_V1 <- lapply(graph_V1, FUN = function(x) data.frame(E(x)$weight)[,1])
  nb_edges <- lapply(edges_V1, length)
  rep_ids <- c()
  for(i in 1:length(idsV1)){
    rep_idi <- lapply(nb_edges[[i]],FUN = function(x) rep(idsV1[[i]],x))
    rep_ids <- c(rep_ids,rep_idi)
  }
  data_edge_V1 <- data.frame(
    sub_ids = factor(unlist(rep_ids)),
    edges = unlist(edges_V1)
  )

  p <-ggplot(data_edge_V1,aes(x = edges,y = sub_ids, fill = ..x..)) +
        geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
        scale_fill_viridis(name = "edge weight",option = 'C') +
        labs(title = paste("Distribution of edge weight",v_id))
        theme_ridges() +
        theme_classic() +
        theme(legend.position = 'none')
  p
}

#' Plot the edges weight distribution for each participant of a group. FA weighted
#'
#'
#' Get fractional anistropy edges weight distribution each participants of a group
#' (V1, V2, V3). Plot it as a ridge plot. Note : You can plot the three ridges side to side
#' by calling this function three times and using gridExtra::grid.arrange.
#'
#' @param dataFA dataframe containing all the FA weighted data
#' @param v_id chr visit id like "V1", "V2"...
#' @export
get_edge_distribution <- function(dataFA,v_id){
  idsV1 <- get_subject_ids(dataFA,v_id)
  graph_V1 <- lapply( idsV1, FUN =  function(x) makeGraph(dataFA,x,v_id,"FA",0))
  edges_V1 <- lapply(graph_V1, FUN = function(x) data.frame(E(x)$weight)[,1])
  nb_edges <- lapply(edges_V1, length)
  rep_ids <- c()
  for(i in 1:length(idsV1)){
    rep_idi <- lapply(nb_edges[[i]],FUN = function(x) rep(idsV1[[i]],x))
    rep_ids <- c(rep_ids,rep_idi)
  }
  data_edge_V1 <- data.frame(
    sub_ids = factor(unlist(rep_ids)),
    edges = unlist(edges_V1)
  )
  data_edge_V1
}



