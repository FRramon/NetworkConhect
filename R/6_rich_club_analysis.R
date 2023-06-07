
#' Compute non normalized rich club coefficient of a graph g
#'
#' @param g a undirect weighted graph
#' @param k the degree threshold in the induced subgraph construction
#' @export
computePhi <- function(g,k){
  wrank <- sort(E(g)$weight,decreasing = T)

  vk <- V(g)[degree(g)>k]
  Gk <- induced_subgraph(g,vids = vk)
  Ek <- length(E(Gk))
  Wk <- sum(E(Gk)$weight)

  swrank <- sum(wrank[1:Ek])

  phi <- Wk/swrank
  phi
  #list("wrank" = wrank,"vk" = vk,"Ek" = Ek,"Wk" = Wk,"phi" = phi)
}

#' Construct a random graph from a graph g
#' Respect the same nb of edges/nodes. And keep the same degree distribution
#' and keep edge weight distribution
#' @param g a undirect weighted graph
#' @export
computeGr <- function(g){
  deg <- degree(g)
  weights <- E(g)$weight
  new_weights <- sample(weights)
  gR <- sample_degseq(deg,method = 'vl') # Conservation of number of nodes & deg distribution
  Random_elist <- data.frame(as_edgelist(g))
  Random_elist$weight <- new_weights
  # print(head(Random_elist))
  gR <- graph_from_data_frame(Random_elist)
  gR
}


#' Compute all phi for all subject in one group (v_id)
#' @param data a dataframe containing all data for analysis
#' @param v_id visit id : "V1","V2","V3"
#' @param metric the metric chosen for the weight in the network. "FBC","FA"..
#' @param tvalue density threshold value. Default .15
#' @export
compute_phi_patients <- function(data,v_id,metric,tvalue,k,nrand){
  subjects <- get_subject_ids(data,v_id)
  P <- c()
  for(s in subjects){
    g <- makeGraph(data,s,v_id,metric,0)
    gt <- sparseThresh(g,tvalue)
    phi <- computePhi(gt,k)

    # gR <- computeGr(g)
    phiR <- c()
    for(nr in 1:nrand){
      gR <- computeGr(gt)
      phir <- computePhi(gR,k)
      #unlist(lapply(X = X, FUN = computePhi,g = gR))
      phiR <- cbind(phiR,phir)
    }
    phirand <- rowMeans(phiR)
    #phir <- computePhi(gR,k)

    phin <- phi/phirand
    P <- c(P,phin)
  }
  print(P)
  print(length(P))
  P

}


#' Do all statistical analysis at each k
#' @param data a dataframe containing all data for analysis
#' @param metric the metric chosen for the weight in the network. "FBC","FA"..
#' @param tvalue density threshold value. Default .15
#' @export
computePhiGroup <- function(dataFBC,metric,t,k,nrand){
  # print(k)
  # print("Loading data...")
  # data_dir <- getDataDir(metric)
  # dataFBC <- read_and_normalize_data(data_dir,metric)
  print("Computing V1...")
  v1 <- compute_phi_patients(dataFBC,"V1",metric,t,k,nrand)
  print("Computing V2...")
  v2 <- compute_phi_patients(dataFBC,"V2",metric,t,k,nrand)
  print("Computing V3...")
  v3 <- compute_phi_patients(dataFBC,"V3",metric,t,k,nrand)

  D <- data.frame(
    richphi = c(v1,v2,v3),
    timepoint = c(rep("V1",length(v1)),rep("V2",length(v2)),rep("V3",length(v3))),
    subject = c(get_subject_ids(dataFBC,"V1"),get_subject_ids(dataFBC,"V2"),get_subject_ids(dataFBC,"V3"))
  )
  print("Fitting LME")
  fit <- lmer(richphi ~ timepoint + (1|subject),data = D)
  list('coeffs' = summary(fit)$coefficients, 'data' = D)
}


#' Main function for the rich club analysis
#' @param data a dataframe containing all data for analysis
#' @param metric the metric chosen for the weight in the network. "FBC","FA"..
#' @param tvalue density threshold value. Default .15
#' @param kmin minimum degree for the range of analysis default 10
#' @param kmax maximum degree for the range of analysis default 50
#' @export
main_rich_club_analysis <- function(data_path,metric = "FBC",tvalue = 0.15 ,kmin = 10,kmax = 50,nrand){

  file_path <- getDataDir(data_path,metric)
  print("Step 1 : Loading data...")
  print(file_path)
  data <-read_and_normalize_data(file_path, metric)
  X <- seq(kmin,kmax,1)
  j=1
  M <- c("V1","V2","V3")
  SD <- c("V1","V2","V3")
  fit <- list()
  for(i in X){
    fit[[j]] <- computePhiGroup(dataFBC = data,metric = metric, t= tvalue, k= i,nrand)
    m1 <- aggregate(x=fit[[j]]$data$richphi, by = list(fit[[j]]$data$timepoint),FUN=mean)
    m2 <- m1$x
    M <- rbind(M,m2)
    sd1 <- aggregate(x=fit[[j]]$data$richphi, by = list(fit[[j]]$data$timepoint),FUN=mean)
    sd2 <- sd1$x
    SD <- rbind(SD,sd2)
    j=j+1
  }

  SD <- c("sdV1","sdV2","sdV3")
  j=1
  for(i in X){
    sd1 <- aggregate(x=fit[[j]]$data$richphi, by = list(fit[[j]]$data$timepoint),FUN=sd)
    sd2 <- sd1$x
    SD <- rbind(SD,sd2)
    j=j+1
  }
  DSD <- as.data.frame(SD)[-1,]
  rownames(DSD) <- X
  colnames(DSD) <- c("sdV1","sdV2","sdV3")
  DSD[] <- lapply(DSD,as.numeric)


  DM <- as.data.frame(M)[-1,]
  rownames(DM) <- X
  colnames(DM) <-  c("mV1","mV2","mV3")
  DM[] <- lapply(DM,as.numeric)

  D <- data.frame(cbind(DM,DSD))

  p <- ggplot(D,aes(x=X)) +
      geom_line(aes(y=mV1,color = "V1")) +
      geom_line(aes(y=mV2,color = "V2")) +
      geom_line(aes(y=mV3,color = "V3")) +
      geom_point(aes(y=mV1,color = "V1")) +
      geom_point(aes(y=mV2,color = "V2")) +
      geom_point(aes(y=mV3,color = "V3")) +
      geom_errorbar(aes(ymin = mV1 - sdV1,ymax = mV1 + sdV1),width = .2) +
      geom_errorbar(aes(ymin = mV2 - sdV2,ymax = mV2 + sdV2),width = .2) +
      geom_errorbar(aes(ymin = mV3 - sdV3,ymax = mV3 + sdV3),width = .2) +

      labs(color = "Curves",x= "degree k",y = "phi normalized") +
      scale_color_manual(values = c("V1" = "lightgreen","V2" = "darkolivegreen3","V3" = "darkgreen")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = c(0.1,0.9),
            legend.background =element_rect(fill = "white",color = "black"),
            legend.key = element_rect(fill = "white"),
            legend.title = element_blank())

  p3 <- list()
  p2 <- list()
  for(j in 1:length(fit)){
    p2[[j]] <- fit[[j]]$coeffs[2,5]
    p3[[j]] <- fit[[j]]$coeffs[3,5]
    print(j)
  }
  list("plot" = p,"pvaluesV1V2" = p2,"pvaluesV1V3" = p3)

}
