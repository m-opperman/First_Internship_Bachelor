###########################################################################################################################################
# Author:   M. Opperman                                                                                                                   #
# Date:     20-11-2018                                                                                                                    #
# Version:  2.0                                                                                                                           #
#                                                                                                                                         #
# Function: Contains the pre proces functions for plotting of heatmaps.                                                                   #
#           These functions are being called from the main_PA.R script                                                                    #                                                                      
###########################################################################################################################################


#' Plot copynumber data in heatmap
#'
#' @description Pre proccesses data to plot the copyNumber in a heatmap and calls heatMap_plot function to visualize copy number data
#'
#' @param data_frame Data frame which contains copy number data, wgd status and primary tumor location per sample
#' @param plotName String which determines the name of the plot and if data should be split up into a smaller group based on wgdStatus
#'
#' @return None
#' @export
#'
#' @examples
#' plotCopyNumber(results_CopyNumber, plotName)
plotCopyNumber <- function(data_frame, plotName){
  uniqueTypes <- unique(data_frame$primaryTumorLocation)
  data_frame  <- wgdNormalization(data_frame, TRUE)
  #data_frame  <- normalization(data_frame, TRUE)

  data_frame1   <- data_frame[which(data_frame$wgdStatus <= 0),]
  data_frameWGD <- data_frame[which(data_frame$wgdStatus >= 1),]

  for(i in c(1:3)){
    if(i == 2){
      data_frame <- data_frame1
      plotName <- "CopyNumber not WGD"
    }else if(i == 3){
      data_frame <- data_frameWGD
      plotName <- "CopyNumber WGD"
    }

    list_data <- lapply(uniqueTypes, function(x){
      subtypeSet <- subset(data_frame, data_frame$primaryTumorLocation == x)
      tmp_list   <- unlist(getMeans(subtypeSet))
    })

    tmp_list <- lapply(uniqueTypes, function(x){
      len <- length(data_frame[which(data_frame$primaryTumorLocation == x),1])
      return(paste0(x, " (", len, ")"))
    })

    heatmap_df <- do.call(cbind, list_data)
    colnames(heatmap_df) <- unlist(tmp_list)
    rownames(heatmap_df) <- getChromArms()

    heatMap_plot(heatmap_df, plotName)
  }
}



#' Plot loss of heterozygosity in heatmap
#'
#' @description Determines per cancer-type the percentage of Loss of heterozygosity for each autosome arm and calls heatmap_plot function to
#' visualize data in a heatmap
#'
#' @param results_MinorAP Data frame which contains minor allele ploidy data, wgd status and primary tumor location per sample
#'
#' @return None
#' @export
#'
#' @examples
#' plotLOH(results_MinorAP)
plotLOH <- function(results_MinorAP){
  uniqueTypes <- unique(results_MinorAP$primaryTumorLocation)
  listLOH <- lapply(uniqueTypes, function(x){
    return(unlist(getPercentage(results_MinorAP[which(results_MinorAP$primaryTumorLocation == x),])))
  })

  tmp_list <- lapply(uniqueTypes, function(x){
    len <- length(results_MinorAP[which(results_MinorAP$primaryTumorLocation == x),1])
    return(paste0(x, " (", len, ")"))
  })

  heatmap_df <- do.call(cbind, listLOH)
  colnames(heatmap_df) <- unlist(tmp_list)
  rownames(heatmap_df) <- getChromArms()
  heatMap_plot(heatmap_df , "loh")
}


#' Plot unsupervised heatmap of copynumber data
#'
#' @description Makes unsupervised heatmap of copy number data, calls plotfunction from heatmapPlots_PA.R
#'
#' @param results_CopyNumber Data frame which contains copy number data, wgd status and primary tumor location per sample
#' @param mut_load_anno Data frame which contains mutational load annotation per sample
#' @param tp character string which defines if data should be split
#'
#' @return None
#' @export
#'
#' @examples
#' preProcesUnsupervisedHeatmap(results_CopyNumber, mut_load_anno, "NonWGD")
#' preProcesUnsupervisedHeatmap(results_CopyNumber, mut_load_anno, "WGD")
#' preProcesUnsupervisedHeatmap(results_CopyNumber, mut_load_anno, "All")
preProcesUnsupervisedHeatmap <- function(results_CopyNumber, mut_load_anno, tp){
  if(tp == "WGD"){
    #data_frame <- wgdNormalization(results_CopyNumber, TRUE)
    data_frame <- wgdNormalization(results_CopyNumber, TRUE)
    data_frame <- data_frame[which(data_frame$wgdStatus >= 1),]
  }else if(tp == "NonWGD"){
    #data_frame <- wgdNormalization(results_CopyNumber, TRUE)
    data_frame <- wgdNormalization(results_CopyNumber, TRUE)
    data_frame <- data_frame[which(data_frame$wgdStatus <= 0),]
  }else if(tp == "All"){
    #data_frame <- wgdNormalization(results_CopyNumber, TRUE)
    data_frame <- wgdNormalization(results_CopyNumber, TRUE)
  }

  annoDF <- as.data.frame(data_frame[,41:42])
  rownames(annoDF) <- rownames(data_frame)

  data_frame <- as.matrix(data_frame[,1:39])
  data_frame <- t(data_frame)

  annoDF     <- merge.data.frame(annoDF, mut_load_anno, by = "row.names")
  typeColors <- c("lightblue","grey","purple","yellow","darkgreen","brown","pink","orange","beige",
                  "cyan","springgreen","black","deepskyblue","salmon","navy","khaki4",
                  "indianred","darkorange4","ivory2","yellowgreen","burlywood3")

  annoColors <- as.data.frame(factor(annoDF$primaryTumorLocation, levels = unique(annoDF$primaryTumorLocation),
                                     labels = typeColors))

  colnames(annoColors)      <- c("typeColor")
  if(tp == "All"){
    annoColors["wgdColor"]  <- factor(annoDF$wgdStatus, levels = unique(annoDF$wgdStatus), labels = c("blue", "grey", "purple", "green"))
  }else{
    annoColors["wgdColor"]  <- factor(annoDF$wgdStatus, levels = unique(annoDF$wgdStatus), labels = c("burlywood3", "black"))
  }
  annoColors["clonalColor"] <- cut(log10(annoDF$mut_load_clonal),right = TRUE, breaks = seq(0,6.1, by = 0.1),
                                   labels = colorpanel(61, low = "white", mid = "yellow", high = "red"))

  heatPlotsUnsupervised(data_frame, annoColors, annoDF, typeColors, tp)
}

