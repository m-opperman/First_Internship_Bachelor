###########################################################################################################################################
# Author:   M. Opperman                                                                                                                   #
# Date:     20-11-2018                                                                                                                    #
# Version:  2.0                                                                                                                           #
#                                                                                                                                         #
# Function: Contains the necessary heatmap plot functions to analyse the ploidy levels of cancer samples.                                 #
#           These functions are being called from the pre_proces_heatmaps_PA.R script                                                     #                                                                      
###########################################################################################################################################


#' Plot function heatmap
#'
#' @description Makes heatmap plot of data with and without a dendogram clustering
#'
#' @param heatmap_df Data frame which contains the data for heatmap plot.
#' @param plotName character string which is used as plot name and determines the the number of breaks for the heatmap
#'
#' @return None
#' @export
#'
#' @examples
#' heatMap_plot(heatmap_df, "CopyNumber not WGD")
#' heatMap_plot(heatmap_df, "CopyNumber WGD")
#' heatMap_plot(heatmap_df, "CopyNumber")
#' heatMap_plot(heatmap_df, "loh")
#' heatMap_plot(heatmap_df, "major")
#' heatMap_plot(heatmap_df, "minor")
#' heatMap_plot(heatmap_df, "major/minor")
heatMap_plot  <- function(heatmap_df, plotName){
  if(plotName == "CopyNumber not WGD" || plotName == "CopyNumber WGD" || plotName == "CopyNumber" || plotName == "loh"){
    if(plotName == "loh"){
      breaks <- seq(from =  0, to = 100, by = 1)
    }else{
      breaks <- seq(from = -1, to =   1, by = 0.05)
    }
    print(heatmap.2(heatmap_df, trace  = "none", density = "none", col  = bluered(length(breaks) -1), dendrogram = "column",
                    breaks = breaks, Rowv = FALSE,  Colv = TRUE,  main = plotName, symkey = FALSE, lhei=c(0.1,1), lwid=c(0.5,1),
                    cexRow = 1.1, cexCol  = 1.5  ,  margins=c(13,10) , cellnote = round(heatmap_df, 2), notecol = "black", notecex = 1))

  }else if(plotName == "major" || plotName == "minor" || plotName == "major/minor"){
    print(ggplot(heatmap_df, aes(x   = plotType, y   = chromArm,  fill = ploidy, group = cancerType)) + geom_tile()+
            scale_fill_gradient2(low = "#0000FF" , mid = "#FFFFFF", high = "#FF0000",na.value = "grey50", midpoint = 1, guide = "colourbar") +
            theme(axis.text.x = element_text(angle = 60, hjust = 1)) + ylab("Autosome Arms") + xlab("Cancer types") +
            ggtitle(plotName) + facet_wrap(~cancerType, scales ="free_x", nrow =2) + geom_text(aes(label = round(ploidy, 1))))
  }
}


#' Plot function unsupervised heatmap
#'
#' @description Plot function for unsupervised heatmap of copy number data
#'
#' @param data_frame Data frame that contains numerical values to plot heatmap
#' @param annoColors Data frame which contains per type of value the right colour code
#' @param annoDF Data frame that contains annotation information
#' @param typeColors Character string that contains 21 color codes (for each cancertype one)
#' @param tp Character string that defines the plot name
#'
#' @return None
#' @export
#'
#' @examples
#' heatPlotUnsupervised(data_frame, annoColors, annoDF, typeColors, "All")
#' heatPlotUnsupervised(data_frame, annoColors, annoDF, typeColors, "Non WGD")
#' heatPlotUnsupervised(data_frame, annoColors, annoDF, typeColors, "WGD")
heatPlotsUnsupervised <- function(data_frame, annoColors, annoDF, typeColors, tp){
  breaks    <- seq(from = -1, to =   1, by = 0.05)
  colnames(annoColors[,1:2]) <- c("WGD Status", "Primary Tumor Location")
  #png(paste0(tp , "_heatmap_type_wgd.png"), width = 3000, height = 2000)
  png(paste0("~/Documents/Packages/ploidyAnalyzer/output/heatmaps/", tp,"_heatmap_type_wgd.png"), width = 1500, height = 750)

  heatmap.plus(data_frame,col  = bluered(length(breaks) -1), scale = "none",
               breaks = breaks, Rowv = NA,  Colv = TRUE, cex =30,  hclustfun = hclust2,
               margins=c(1,20), cexRow = 1.5,    cexCol  = 0.1,  keep.dendro = FALSE,
               useRaster=TRUE, ColSideColor = as.matrix(annoColors[,1:2]))

  if(tp == "All"){
    legend("bottomright", 5,5, legend = c(-1,0,1,2), fill = c("red","blue", "grey", "purple"), cex=0.85)
  }else{
    legend("bottomright", 5,5, legend = unique(annoDF$wgdStatus),fill=c("burlywood3", "black"),cex=0.85)
  }
  legend("topright",5,5, legend= unique(annoDF$primaryTumorLocation), fill = typeColors,cex=0.85)
  legend("right",   5,5, legend= c("-1: Loss","  0: No gain or loss","  1: Gain"), fill = c("blue", "white", "red"), cex = 0.85)

  dev.off()
}


#' Plot function for unsupervised loss of heterozygosity heatmap
#'
#' @description Plot function for unsupervised heatmap of results minor alle ploidy data
#'
#' @param results_MinorAP Data frame which contains minor allele ploidy data, wgd status and primary tumor location per sample
#' @param mut_load_anno Data frame which contains mutational load annotations per sample
#'
#' @return None
#' @export
#'
#' @examples
#' lohUnsupervisedHeatmap(results_MinorAP, mut_load_anno)
lohUnsupervisedHeatmap <- function(results_MinorAP, mut_load_anno){
  breaks     <- seq(from = 0, to = 2, by = 0.05)
  data_frame <- results_MinorAP[,1:39]
  data_frame[data_frame > 2] <- 2
  data_frame <- t(data_frame)
  annoDF     <- results_MinorAP[,41:42]
  rownames(annoDF) <- rownames(results_MinorAP)

  typeColors <- c("lightblue","grey","purple","yellow","darkgreen","brown","pink","orange","beige",
                  "cyan","springgreen","black","deepskyblue","salmon","navy","khaki4",
                  "indianred","darkorange4","ivory2","yellowgreen","burlywood3")

  annoColors <- as.data.frame(factor(annoDF$primaryTumorLocation, levels = unique(results_MinorAP$primaryTumorLocation),
                                     labels = typeColors))

  colnames(annoColors)   <- c("typeColor")
  annoColors["wgdColor"] <- factor(annoDF$wgdStatus, levels = c(-1,0,1,2), labels = c("red","blue", "grey", "purple"))
  #png("loh_unsupervised.png", width = 3000, height = 2000)
  png("~/Documents/Packages/ploidyAnalyzer/output/heatmaps/loh_unsupervised.png", width = 1500, height = 750)

  heatmap.plus(as.matrix(data_frame),col = bluered(length(breaks) -1), scale = "none", cex.main = 10,
               main = "Unsupervised Clustering Of Normalized Minor Allele Ploidy Data",
               Rowv = NA,  Colv = TRUE, cex =3, hclustfun = hclust2, breaks = breaks,
               margins=c(1,20), cexRow = 1.5,    cexCol  = 0.1,  keep.dendro = FALSE,
               useRaster=TRUE, ColSideColor = as.matrix(annoColors))

  legend("right",      5,5,  legend = c("0: LOH","1: No gain or loss","2: Gain"), fill = c("blue", "white", "red"), cex=0.85)
  legend("bottomright",5,5,legend = c("-1: Haploid"," 0: Non WGD"," 1: One WGD"," 2: Two WGDs"),fill=c("red","blue","grey","purple"),cex=0.85)
  legend("topright",   5,5,legend = unique(annoDF$primaryTumorLocation),fill= typeColors,cex=0.85)

  dev.off()

}






