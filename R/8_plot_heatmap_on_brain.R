

#' @export
get_data_visit <- function(metric,eval,v_id, thresh_method,tvalue){
#  csv_res <- read.csv(paste("/home/imabrain/Documents/Brain_networks/vAllWeights/Results/",metric,"/",v_id,"_nodal_",eval,wantedDensity,".csv",sep=""))
  print(paste("/home/imabrain/Documents/Brain_networks/vAllWeights/Results/",metric,"/",v_id,"_nodal_",eval,thresh_method,as.character(tvalue),".csv",sep=""))
  csv_res <- read.csv(paste("/home/imabrain/Documents/Brain_networks/vAllWeights/Results/",metric,"/",v_id,"_nodal_",eval,thresh_method,as.character(tvalue),".csv",sep=""))

  c1 <- rowMeans(csv_res[,-1])
  res1 <- data.frame(cbind(csv_res[,1],c1))
  colnames(res1) <- c('No','value')

  lut <- getLUT()

  res_labelname <- join(res1,lut,by = 'No')
  res_labelname$Y <- v_id
  res <- res_labelname[,c(3,4,2)]
  colnames(res) <- c('X','Y','Z')
  res
}

#' @export
get_res <- function(metric,eval, thresh_method,tvalue){
  res1 <- get_data_visit(metric,eval,"V1",thresh_method,tvalue)
  res2 <- get_data_visit(metric,eval,"V2",thresh_method,tvalue)
  res3 <- get_data_visit(metric,eval,"V3",thresh_method,tvalue)

  res <- as.data.frame(rbind(res1,res2,res3))
}

#' @export
get_minmax <- function(res){
  minV <- min(res$Z)
  maxV <- max(res$Z)
  list("min" = minV, "max" = maxV)
}

#' @export
main_cortical_values <- function(metric,eval,thresh_method,tvalue){

  res <- get_res(metric,eval,thresh_method,tvalue)


   # return(lims)
  # csv res du type comp
  #Ici on split le df comp en zone cort et subcort, l/r
  res_cort <- split_zone_hemi(res)

  min_max_cortl <- get_minmax(res_cort$cortl)
  min_max_cortr <- get_minmax(res_cort$cortr)
  min_cort <- min(min_max_cortl$min,min_max_cortr$min)
  max_cort <- max(min_max_cortl$max,min_max_cortr$max)
  # Ici on récupère les noms des labels utilisés dans ggseg
  cort_atlas_labels <- get_cortical_labels()
  # On opère la jointure entre les deux :
  cortical_res <- join_res_labels(res_cort,cort_atlas_labels)
  #return(cortical_res)
  # Les resultats de la jointure
  res_cortl <- cortical_res$cortl
  res_cortr <- cortical_res$cortr
  # Definition des data à plot avec tibble
  data_cortl <- tibble(
    region = as.character(res_cortl$region),
    centrality = res_cortl$Z,
    Time = as.character(res_cortl$Y)
  )

  data_cortr <- tibble(
    region = as.character(res_cortr$region),
    centrality = res_cortr$Z,
    Time = as.character(res_cortr$Y)
  )

  #return(data_cortr)
  # definition des figures
  figCortl <- data_cortl %>%
    group_by(Time) %>%
    ggseg(atlas = "desterieux",hemi = "left",position = "stacked",
          mapping = aes(fill = centrality),
          show.legend = TRUE, colour = "black") +
    facet_wrap(~Time, ncol= 3) +
    scale_fill_gradientn(colours = rev(rainbow(20*10,start = 0/6, end = 4/6)),
                         limits = c(min_cort,max_cort)) +
    labs(title = "Left Hemisphere")

  figCortr <- data_cortr %>%
    group_by(Time) %>%
    ggseg(atlas = "desterieux",hemi = "right",position = "stacked",
          mapping = aes(fill = centrality),
          show.legend = TRUE, colour = "black") +
    facet_wrap(~Time, ncol= 3) +
    scale_fill_gradientn(colours = rev(rainbow(20*10,start = 0/6, end = 4/6)),
                         limits = c(min_cort,max_cort)) +
    labs(title = "Right Hemisphere")

  # plot une ligne deux colonnes
  #gridExtra::grid.arrange(figCortl,figCortr,ncol = 1,nrow =2 )
  list("cortl" = figCortl, "cortr" = figCortr)

}

#' @export
main_subcortical_values <- function(metric,eval,thresh_method,tvalue){

  res <- get_res(metric,eval,thresh_method,tvalue)
  # csv res du type comp
  #Ici on split le df comp en zone cort et subcort, l/r
  res_cort <- split_zone_hemi(res)
  min_max_subl <- get_minmax(res_cort$subl)
  min_max_subr <- get_minmax(res_cort$subr)
  min_sub <- min(min_max_subl$min,min_max_subr$min)
  max_sub <- max(min_max_subl$max,min_max_subr$max)
  print(max_sub)
  print(min_sub)
  #return(min_sub)
  # Ici on récupère les noms des labels utilisés dans ggseg
  sub_atlas_labels <- get_subcortical_labels()
  # On opère la jointure entre les deux :
  sub_res <- join_res_sub_labels(res_cort, sub_atlas_labels)
  # Les resultats de la jointure
  res_subl <- sub_res$subl
  res_subr <- sub_res$subr
  # Definition des data à plot avec tibble
  data_subl <- tibble(
    region = as.character(res_subl$region),
    centrality = res_subl$Z,
    Time = as.character(res_subl$Y)
  )

  data_subr <- tibble(
    region = as.character(res_subr$region),
    centrality = res_subr$Z,
    Time = as.character(res_subr$Y)
  )

  # definition des figures
  figSubl <- data_subl %>%
    group_by(Time) %>%
    ggseg(atlas = "aseg",hemi = "left",position = "stacked",
          mapping = aes(fill = centrality),
          show.legend = TRUE, colour = "black") +
    facet_wrap(~Time, ncol= 3) +
    scale_fill_gradientn(colours = rev(rainbow(20*10,start = 0/6, end = 4/6)),
                         limits = c(min_sub,max_sub)) +
    labs(title = "Left Hemisphere")

  figSubr <- data_subr %>%
    group_by(Time) %>%
    ggseg(atlas = "aseg",hemi = "right",position = "stacked",
          mapping = aes(fill = centrality),
          show.legend = TRUE, colour = "black") +
    facet_wrap(~Time, ncol= 3) +
    scale_fill_gradientn(colours = rev(rainbow(20*10,start = 0/6, end = 4/6)),
                         limits = c(min_sub,max_sub)) +
    labs(title = "Right Hemisphere")

  # plot une ligne deux colonnes
  #gridExtra::grid.arrange(figSubl,figSubr,ncol = 1,nrow =2 )
  list("subl" = figSubl,"subr" = figSubr)
}

#' @export
plot_values_on_brain_map <- function(metric,eval,thresh_method,tvalue){
  templates <- get_templates()
  res_cort <- main_cortical_values(metric,eval,thresh_method,tvalue)
  res_subcort <- main_subcortical_values(metric,eval,thresh_method,tvalue)
  pagePlot <- gridExtra::grid.arrange(res_cort$cortl,res_cort$cortr,res_subcort$subl,res_subcort$subr,
                          ncol = 2,
                          nrow = 2,
                          top = grid::textGrob(paste("Nodal",eval,"in",metric,"weighted network - ",as.character(tvalue),thresh_method),gp = grid::gpar(fontsize = 20, font =3)))
  print(pagePlot)
  ggsave(filename = paste("/home/imabrain/Documents/Brain_networks/vAllWeights/Results/",metric,"/Local/nodal_",eval,tvalue,".png",sep=""),pagePlot)
}
