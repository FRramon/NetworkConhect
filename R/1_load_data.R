# Author Francois Ramon

###############################################################################################
# FIND DATA DIRECTORY, LOAD DATA GIVEN THE METRIC CHOICE, FORMAT : FROM,TO,WEIGHT
###############################################################################################

#' Find the path containing the .xlsx file of the asked weighting scheme.
#'
#' @param choice a chr for the weighting scheme : "FA", "FBC", "GFA"...
#' @returns chr the path containing the .xlsx file
#' @export
#'
#'

getDataDir<-function(data_path, choice = c('FBC','ODI','GFA','FA','Fintra','PearsonCorrel')){

  choice<-match.arg(choice)
 # data_path <- '/home/imabrain/Documents/GraphConhect/data'
  files <- list.files(data_path,recursive=FALSE)
  filename = paste('stats_diffusion_metric_in_fullWM_',choice, '.xlsx',sep = "")
  print(filename)
  if(filename %in% files){
    paste(data_path,'/',filename,sep="")
   # print('ok!')
  }else{
    print("file not found")
  }
}

#' Create a dataframe with columns "subject_id","visit_id","from","to","weight" for the chosen weighting scheme
#'
#' @param data_path chr indicating the path of the "raw" data
#' @param metric a chr for the weighting scheme : "FA", "FBC", "GFA"...
#' @returns the newly created dataframe
#' @export
read_and_normalize_data <- function(data_path,WM_metric){
  data_full <- read_excel(data_path)
  if(WM_metric=='FBC'){
    # drois <- separateRois(data_full)
    # weight: 9 for fiber density, 8 for distance, 5 for streamline count
    #print("select streamline count")
    df <- as.data.frame(cbind(data_full[,1],data_full[,2],data_full[,6],data_full[,7],data_full[,5]))
    colnames(df) <- c('subject_id','visit_id','from','to','weight')
    df
  }else{
    df <- as.data.frame(cbind(data_full[,1],data_full[,2],data_full[,12],data_full[,13],data_full[,10]))
    colnames(df) <- c('subject_id','visit_id','from','to','weight')
    df
  }
}

#
# read_and_normalize_data <- function(data_path,WM_metric=c('FBC','GFA','FA','ODI','Fintra','FA')){
#   data_full <- read_excel(data_path)
#   if(WM_metric=='FBC'){
#     drois <- separateRois(data_full)
#     df <- as.data.frame(cbind(drois[,1],drois[,2],drois[,6],drois[,7],drois[,5]))
#     colnames(df) <- c('subject_id','visit_id','from','to','weight')
#     df
#   }else{
#     df <- as.data.frame(cbind(data_full[,1],data_full[,2],data_full[,12],data_full[,13],data_full[,10]))
#     colnames(df) <- c('subject_id','visit_id','from','to','weight')
#     df
#   }
# }


#' Create a dataframe from freesurfer lut.
#'
#' @returns dataframe of two columns : No for the region id, labelname for the region label
#'
#'
#' @export
getLUT <- function(data_path){
  #txt_path <- paste(data_path,"freesurfer_lut.txt",sep="")
  doc <- read.table(data_path)
  colnames(doc) <- c('No','labelname','r','g','b','A')
  df<-as.data.frame(doc)
  df[,c('No','labelname')]
}

#' Get the list of subject participating visit i
#'
#' @param data dataframe (created by 'read_and_normalize')
#' @param v_id chr visit id, like "V1", "V2", "V3".
#' @returns list
#' @export
get_subject_ids <- function(data,v_id){
  df <- subset(data,visit_id == v_id)
  unique(df[,1])
}


