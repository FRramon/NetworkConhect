
# Enregitrer ces données dans un csv file, et les load ici (donc faire le write.csv dans 10)

#' Split zones in the brain : left/rigth hemisphere, and cortical/subcortical
#' @param csv_res a dataframe with columns X, Y , Z, coming from local analysis
#' @returns list of sub dataframe, for each zone
#' @export
split_zone_hemi <- function(csv_res){

  # parcours lignes label : si ctx alors cortex, sinon sub
  labels <- unique(csv_res$X)
  cortical_areas <- grep("ctx",labels,value  = TRUE)
  subcortical_areas <- grep("Left|Right",labels, value = TRUE)
  cortical_left <- data.frame(X = grep("lh",cortical_areas,value = TRUE))
  cortical_right <- data.frame(X = grep("rh",cortical_areas, value = TRUE))

  subcortical_left <- data.frame(X = grep("Left", subcortical_areas, value = TRUE))
  subcortical_right <- data.frame(X = grep("Right", subcortical_areas, value = TRUE))
  print(subcortical_right)
  res_cortl <- join(cortical_left,csv_res, by = "X")
  res_cortr <- join(cortical_right,csv_res, by = "X")
  res_subl <- join(subcortical_left,csv_res, by = "X")
  res_subr <- join(subcortical_right,csv_res, by = "X")

  list("cortl" = res_cortl,"cortr" = res_cortr,"subl" = res_subl, "subr" = res_subr)
}

#' Table of region names equivalence in the cortical region
#'
#' Create a equivalence table between region names of destrieux atlas in the package ggseg, and region names used in freesurfer.For cortical regions
#' @returns list of two tables for left and right hemisphere. tables of equivalence
#' @export
get_cortical_labels <- function(){
  destr<- ggseg (atlas = "desterieux",
                 mapping = aes(fill = label)) +
    theme(legend.justification = c(1, 0),
          legend.position = "bottom" ,
          legend.text = element_text(size = 5)) +
    guides(fill = guide_legend(ncol = 3))

  desterieux_labels_lh <- lapply(unique(destr$data$label)[1:75], FUN = function(x) paste("ctx_",x,sep = ""))
  L_equivalent <- cbind(unique(destr$data$region),desterieux_labels_lh)
  L_equivalent <- data.frame(L_equivalent)
  colnames(L_equivalent) <- c("region","X")

  desterieux_labels_rh <- lapply(unique(destr$data$label)[1:75], FUN = function(x) paste("ctx_r",substring(x,2),sep = ""))
  R_equivalent <- cbind(unique(destr$data$region),desterieux_labels_rh)
  R_equivalent <- data.frame(R_equivalent)
  colnames(R_equivalent) <- c("region","X")

  list("Left" = L_equivalent, "Right" = R_equivalent)
}


#' Table of region names equivalence in the subcortical region
#'
#' Create a equivalence table between region names of destrieux atlas in the package ggseg, and region names used in freesurfer.For subcortical regions
#' @returns list of two tables for left and right hemisphere. tables of equivalence
#' @export
get_subcortical_labels <- function(){
  subatlas <- ggseg (atlas = "aseg",
                     mapping = aes(fill = region))

  sublabel_r <- grep("Right",unique(subatlas$data$label),value = TRUE)
  subregion_r <- c("thalamus proper","lateral ventricle","putamen","amygdala","ventral DC","hippocampus","pallidum","caudate","cerebellum cortex","cerebellum white matter")
  sublabel_r[1] <- paste(sublabel_r[1],"*",sep = "")

  R_equivalent <- cbind(subregion_r,sublabel_r)
  R_equivalent <- data.frame(R_equivalent)
  colnames(R_equivalent) <- c("region","X")

  sublabel_l <- grep("Left",unique(subatlas$data$label),value = TRUE)
  subregion_l <- c("thalamus proper","hippocampus","lateral ventricle","putamen","ventral DC","amygdala","pallidum","caudate")
  sublabel_l[1] <- paste(sublabel_l[1],"*",sep = "")

  L_equivalent <- cbind(subregion_l,sublabel_l)
  L_equivalent <- data.frame(L_equivalent)
  colnames(L_equivalent) <- c("region","X")

  list("Left" = L_equivalent,"Right" = R_equivalent)

}


#' Joint of results and region names. cortical regions
#'
#' Join the two tables by the freesurfer name (unique key), because ggseg cannot plot with freesurfer names
#' @param rest_cort dataframe created in split_zone_hemi
#' @param atlas_subcort dataframe of the tables of equivalence for cortical regions
#' @returns list of two tables for left and right hemisphere. tables of equivalence
#' @export
join_res_labels <- function(res_cort,atlas_cort){
  res_cortl <- join(res_cort$cortl,atlas_cort$Left,by = "X")
  res_cortr <- join(res_cort$cortr,atlas_cort$Right,by = "X")
  list("cortl" = res_cortl, "cortr" = res_cortr)
}


#' Joint of results and region names. subcortical regions
#'
#' Join the two tables by the freesurfer name (unique key), because ggseg cannot plot with freesurfer names
#' @param rest_cort dataframe created in split_zone_hemi
#' @param atlas_subcort dataframe of the tables of equivalence for subcortical regions
#' @returns list of two tables for left and right hemisphere. tables of equivalence
#' @export
join_res_sub_labels <- function(res_cort,atlas_subcort){

  res_subl <- inner_join(res_cort$subl,atlas_subcort$Left,by = "X")
  res_subr <- inner_join(res_cort$subr,atlas_subcort$Right,by = "X")
  list("subl" = res_subl, "subr" = res_subr)
}


#' Plot Kruskal Wallis result on cortical regions
#'
#' Use split_zone_hemi, get_(sub)cortical_labels, to get the values to plot. Use ggseg to plot the values on the destrieux atlas.
#' Plot only V1 -V2 and V1 - V3
#' @param metric chr the weighting scheme such as "FA", "GFA"
#' @param eval chr the global network topology measure : like clust_coeff, global_eff
#' @param thresh_method chr the method of threshoding. Whether "threshol" or "density"
#' @param threshold float the threshold/density to be applied on each network in the analysis.
#' @returns two figures : cortical and subcortical figures
#' @export
main_cortical_p <- function(metric,eval,thresh_method,tvalue){

  csv_res <- read.csv(paste("/home/imabrain/Documents/Brain_networks/vAllWeights/Results/",metric,"/nodal_",eval,thresh_method,as.character(tvalue),".csv",sep=""))
  # csv res du type comp
  #Ici on split le df comp en zone cort et subcort, l/r
  res_cort <- split_zone_hemi(csv_res)
  # Ici on récupère les noms des labels utilisés dans ggseg
  cort_atlas_labels <- get_cortical_labels()
  # On opère la jointure entre les deux :
  cortical_res <- join_res_labels(res_cort,cort_atlas_labels)
  # Les resultats de la jointure
  res_cortl <- cortical_res$cortl
  res_cortr <- cortical_res$cortr
  # Definition des data à plot avec tibble
  data_cortl <- tibble(
    region = as.character(res_cortl$region),
    p = res_cortl$Z,
    Time = as.character(res_cortl$Y)
  )

  data_cortr <- tibble(
    region = as.character(res_cortr$region),
    p = res_cortr$Z,
    Time = as.character(res_cortr$Y)
  )

  # definition des figures
  figCortl <- data_cortl %>%
    group_by(Time) %>%
    ggseg(atlas = "desterieux",hemi = "left",position = "stacked",
          mapping = aes(fill = p),
          show.legend = TRUE, colour = "black") +
    facet_wrap(~Time, ncol= 3) +
    scale_fill_gradientn(colours = rainbow(20*10,start = 0/6, end = 4/6),
                         limits = c(-1,1)) +

    labs(title = "Left Hemisphere")

  figCortr <- data_cortr %>%
    group_by(Time) %>%
    ggseg(atlas = "desterieux",hemi = "right",position = "stacked",
          mapping = aes(fill = p),
          show.legend = TRUE, colour = "black") +
    facet_wrap(~Time, ncol= 3) +
    scale_fill_gradientn(colours = rainbow(20*10,start = 0/6, end = 4/6),
                         limits = c(-1,1)) +

    labs(title = "Right Hemisphere")

  # plot une ligne deux colonnes
  #gridExtra::grid.arrange(figCortl,figCortr,ncol = 1,nrow =2 )
  list("cortl" = figCortl, "cortr" = figCortr)

}

#' @export
main_subcortical_p <- function(metric,eval,thresh_method,tvalue){

  csv_res <- read.csv(paste("/home/imabrain/Documents/Brain_networks/vAllWeights/Results/",metric,"/nodal_",eval,thresh_method,as.character(tvalue),".csv",sep=""))
  # csv res du type comp
  #Ici on split le df comp en zone cort et subcort, l/r
  res_cort <- split_zone_hemi(csv_res)
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
    p = res_subl$Z,
    Time = as.character(res_subl$Y)
  )

  data_subr <- tibble(
    region = as.character(res_subr$region),
    p = res_subr$Z,
    Time = as.character(res_subr$Y)
  )

  # definition des figures
  figSubl <- data_subl %>%
    group_by(Time) %>%
    ggseg(atlas = "aseg",hemi = "left",position = "stacked",
          mapping = aes(fill = p),
          show.legend = TRUE, colour = "black") +
    facet_wrap(~Time, ncol= 3) +
    scale_fill_gradientn(colours = rainbow(20*10,start = 0/6, end = 4/6),
                         limits = c(-1,1)) +

    labs(title = "Left Hemisphere")

  figSubr <- data_subr %>%
    group_by(Time) %>%
    ggseg(atlas = "aseg",hemi = "right",position = "stacked",
          mapping = aes(fill = p),
          show.legend = TRUE, colour = "black") +
    facet_wrap(~Time, ncol= 3) +
    scale_fill_gradientn(colours = rainbow(20*10,start = 0/6, end = 4/6),
                         limits = c(-1,1)) +
    labs(title = "Right Hemisphere")

  # plot une ligne deux colonnes
  #gridExtra::grid.arrange(figSubl,figSubr,ncol = 1,nrow =2 )
  list("subl" = figSubl,"subr" = figSubr)
}

#' @export
get_templates <- function(){
  cort <- ggseg(atlas = "desterieux", mapping = aes(fill = region),colour = "black") +
    theme(#legend.justification = c(1,0),
      legend.position = "bottom",
      legend.text = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4))


  subcort <- ggseg(atlas = "aseg", mapping = aes(fill = region),colour = "black") +
    theme(#legend.justification = c(1,0),
      legend.position = "bottom",
      legend.text = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4))

  gridExtra::grid.arrange(cort,subcort,ncol = 2, nrow =1)
  list("cort" = cort, "subcort" = subcort)
}


# MAIN

#' plot significant values on brain map
#' @export
plot_p_on_brain_map <- function(metric,eval,thresh_method,tvalue){
  res_cort <- main_cortical_p(metric,eval,thresh_method,tvalue)
  res_subcort <- main_subcortical_p(metric,eval,thresh_method,tvalue)
  pagePlot <- gridExtra::grid.arrange(res_cort$cortl,res_cort$cortr,res_subcort$subl,res_subcort$subr,
                          ncol = 2,
                          nrow = 2,
                          top = grid::textGrob(paste("Nodal",eval,"in",metric,"weighted network - ",tvalue,thresh_method),gp = grid::gpar(fontsize = 20, font =3)))
  print(pagePlot)
  ggsave(filename = paste("/home/imabrain/Documents/Brain_networks/vAllWeights/Results/",metric,"/Local/nodal_",eval,tvalue,".png",sep=""),pagePlot)

}
