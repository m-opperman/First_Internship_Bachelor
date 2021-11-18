###########################################################################################################################################
# Author:   M. Opperman                                                                                                                   #
# Date:     20-11-2018                                                                                                                    #
# Version:  2.0                                                                                                                           #
#                                                                                                                                         #
# Function: Contains the necessary bar / box plot functions to analyse the ploidy levels of cancer samples.                               #
#           These functions are being called from the main_PA.R script                                                                    #                                                                      
###########################################################################################################################################


#' Duplicated autosome arms plot
#'
#' @description Plots number of autosome arms duplicated on xaxis and the number of samples with the same number on yaxis.
#'
#' @param results_MajorAP Data frame with major allele ploidy data for each autosome arm of each sample and annotation information
#'
#' @return None
#' @export
#'
#' @examples
#' duplicated_arms_plot(results_MajorAP)
duplicated_arms_plot <- function(results_MajorAP){
  output <- lapply(1:length(results_MajorAP$`1p`), function(i){
    ones <- sum(results_MajorAP[i,1:39] == 1)
    zero <- sum(results_MajorAP[i,1:39] == 0)
    rest <- 39 - (ones + zero)
    return(rest)
  })
  x_vector <- (unlist(output))
  y_vector <- rownames(results_MajorAP)
  plot_df  <- as.data.frame(cbind(x_vector, y_vector))

  print(ggplot(data= plot_df, aes(x= factor(x_vector, levels = c(0:39)), fill = factor(results_MajorAP$wgdStatus, levels= c(-1:2))))+
          geom_bar(position = "stack", width = 0.9, stat = "count") +
          ggtitle("Duplicated Autosome Arms") +
          xlab("Number Of Duplicated Autosome Arms")    + ylab("Number Of Samples") +
          guides(fill=guide_legend(title="WGD Status")) +
          scale_fill_manual(labels = c("-1: Haploid", " 0: Non WGD", " 1: One WGD", " 2: Two WGDs" ),values=c("#56B4E9","#8956E9","#999999", "#E69F00")))
}


#' Outlier barplot
#'
#' @description Determines for each high value (outlier) the autosome arm and cancertypes and plots it in barplot
#'
#' @param data_frame Can be either one of the results_ data frames
#' @param treshold Treshold of when a value is considered as outlier
#'
#' @return None
#' @export
#'
#' @examples
#' outliers(results_CopyNumber, 10)
#' outliers(results_MajorAP,     8)
#' outliers(results_MinorAP,     4)
outliers <- function(data_frame, treshold){
  chromArms <- getChromArms()
  samples <- c()
  counter <- 0
  for(x in 1:39){
    if(sum(data_frame[,x] > treshold) != 0){
      tmpSample <- rownames(data_frame[which(data_frame[,x] > treshold),])
      samples <- c(samples, tmpSample)

      if(counter == 0){
        plotDF <- data.frame(subset(data_frame, rownames(data_frame) %in% tmpSample))
        plotDF["chrom"] <- rep(chromArms[x], each = length(tmpSample))
      }else{
        tmp <- data.frame(subset(data_frame, rownames(data_frame) %in% tmpSample))
        tmp["chrom"] <- rep(chromArms[x], each = length(tmpSample))
        plotDF <- rbind(plotDF, tmp)
      }
      counter <- counter + 1
    }
  }
  print(ggplot(data= plotDF, aes(x= factor(chrom, levels = chromArms), fill = primaryTumorLocation))+
          geom_bar(position = "stack", width = 0.9, stat = "count") +
          guides(fill=guide_legend(title="Primary Tumor Location")) +
          ggtitle("Distribution Of Copynumber Outliers") + labs(subtitle = paste0("Outlier if copynumber > ", treshold)) +
          xlab("Autosome Arms")    + ylab("Number Of Outliers"))
}



#' Loss of Heterozygosity distribution barplot
#'
#' @description Determines LOH distribution among autosomearms and cancer types an plots it in a barplot
#'
#' @param data_frame data frame with minor allele ploidy data per autosome arm, per sample
#'
#' @return None
#' @export
#'
#' @examples
#' lohDistribution(results_MinorAP)
lohDistribution <- function(data_frame){
  chromArms <- getChromArms()
  samples <- c()
  counter <- 0
  for(x in 1:39){
    if(sum(data_frame[,x] == 0) != 0){
      tmpSample <- rownames(data_frame[which(data_frame[,x] == 0),])
      samples <- c(samples, tmpSample)
      if(counter == 0){
        plotDF <- data.frame(subset(data_frame, rownames(data_frame) %in% tmpSample))
        plotDF["chrom"] <- rep(chromArms[x], each = length(tmpSample))
      }else{
        tmp <- data.frame(subset(data_frame, rownames(data_frame) %in% tmpSample))
        tmp["chrom"] <- rep(chromArms[x], each = length(tmpSample))
        plotDF <- rbind(plotDF, tmp)
      }
      counter <- counter + 1
    }
  }
  print(ggplot(data = plotDF, aes(x= factor(chrom, levels = chromArms), fill = primaryTumorLocation))+
          geom_bar(position = "stack", width = 0.9, stat = "count") +
          geom_hline(yintercept = 694, linetype = "dashed", color = "black", size = 1) +
          ggtitle("LOH Distribution") + labs(subtitle = "All Samples")+
          xlab("Autosome Arms")    + ylab("Number Of LOH"))
}


#' Whole genome duplication distribution per cancer type barplot
#'
#' @description For each cancer type the contribution of whole genome duplication statuses is determined and visualized
#'
#' @param results_CopyNumber Data frame which contains copy number data for each autosome per sample
#'
#' @return None
#' @export
#'
#' @examples
#' cancerTypes(results_CopyNumber)
cancerTypes <- function(results_CopyNumber){
  options(digits = 2)
  uniqueTypes <- unique(results_CopyNumber$primaryTumorLocation)
  return_list <- lapply(uniqueTypes, function(x){
    total  <- length(results_CopyNumber[which(results_CopyNumber$primaryTumorLocation == x),1])
    tmpSub <- results_CopyNumber[which(results_CopyNumber$primaryTumorLocation == x),]
    wgd    <- length(tmpSub[which(tmpSub$wgdStatus >= 1),1])
    return(c(x, total, wgd))
  })

  wgdFreqDF <- as.data.frame(do.call(rbind,return_list))
  colnames(wgdFreqDF) <- c("primaryTumorLocation", "total_samples", "number_of_wgd")
  wgdFreqDF["percentage"] <- as.numeric(wgdFreqDF$number_of_wgd) / as.numeric(wgdFreqDF$total_samples)
  overalPercentage <- mean(wgdFreqDF$percentage)

  print(ggplot(data = results_CopyNumber, aes(x= primaryTumorLocation, fill = factor(wgdStatus, levels = c(-1:2)))) +
          geom_bar(position = "fill", width = 0.9, stat = "count") +
          ggtitle("Whole Genome Duplication Status Distribution Per Primary Tumor Location") + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
          guides(fill=guide_legend(title="WGD Status")) + xlab("Primary Tumor Location") + ylab("") +
          scale_fill_manual(labels = c("-1: Haploid", " 0: Non WGD", " 1: One WGD", " 2: Two WGDs" ),values=c("#56B4E9","#999999","#8956E9", "#E69F00")) +
          scale_y_continuous(breaks = sort(c(seq(0, 1, length.out=5), overalPercentage))) +
          geom_hline(yintercept = overalPercentage, linetype = "dashed", color = "red", size = 1))
}


#' Distribution of gains, losses and aneuploidyScore per WGD status (Violin plot)
#'
#' @description Plots distribution of gains and losses per wgd Status and primaryTumorLocation
#'
#' @param aneuploidyScoreDF_wgd Data frame which contains gains, losses and primary wgdStatus for each sample.
#'
#' @return None
#' @export
#'
#' @examples
#' violinAneuploidy(aneuploidyScoreDF_wgd)
#' violinAneuploidy(aneuploidyScoreDF)
violinAneuploidy <- function(aneuploidyScoreDF_wgd){
  View(aneuploidyScoreDF_wgd)
  print(ggplot(aneuploidyScoreDF_wgd, aes(x = factor(wgdStatus, levels= c(-1:2)), y = aneuploidyScore  )) +
          geom_violin(scale = "count", draw_quantiles = c(0.25,0.5,0.75),  position=position_dodge(width=0.5)) +
          ggtitle("Aneuploidy Score Distribution Per Whole Genome Duplication Status") +
          xlab("Whole Genome Duplication Status") + ylab("Aneuploidy Score") + labs(subtitle = paste0("Normalization: WGD")))

  aneuploidyScoreDF <- aneuploidyScoreDF_wgd[,1:2]
  aneuploidyScoreDF <- cbind(aneuploidyScoreDF, aneuploidyScoreDF_wgd[,4])
  colnames(aneuploidyScoreDF) <- c("gains", "losses", "wgdStatus")

  plotDF <- melt(aneuploidyScoreDF, id.vars = "wgdStatus")
  plotDF$wgdStatus[which(plotDF$wgdStatus == -1)] <- "-1: Haploid"
  plotDF$wgdStatus[which(plotDF$wgdStatus ==  0)] <- " 0: Non WGD"
  plotDF$wgdStatus[which(plotDF$wgdStatus ==  1)] <- " 1: One WGD"
  plotDF$wgdStatus[which(plotDF$wgdStatus ==  2)] <- " 2: Two WGDs"

  print(ggplot(plotDF, aes(x  = as.factor(variable), y = value, group = variable)) + geom_violin()+
          xlab("Whole Genome Duplication Status") + ggtitle("Distribution Of Gains And Losses Per Whole Genome Duplication Status") +
          ylab("Number Of Gains or Losses") + facet_wrap(~wgdStatus, scales ="free_x"))
}


























