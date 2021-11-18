###########################################################################################################################################
# Author:   M. Opperman                                                                                                                   #
# Date:     20-11-2018                                                                                                                    #
# Version:  2.0                                                                                                                           #
#                                                                                                                                         #
# Function: Contains the necessary support functions to analyse the ploidy levels of cancer samples.                                      #
#           These functions are being called from the main_PA.R script                                                                    #                                                                      #
###########################################################################################################################################


#' Add WGD state to data frame
#'
#' @description Determines the Whole Genome Duplication (WGD) status for each sample and adds it to data_frame
#'
#' @param results_MajorAP Data frame with major allele ploidy data for each autosome arm of each sample and annotation information
#' @param results_MinorAP Data frame with minor allele ploidy data for each autosome arm of each sample and annotation information
#'
#' @return results_MajorAP
#' @export
#'
#' @examples
#' addSampleState(results_MajorAP, results_MinorAP)
addSampleState <- function(results_MajorAP, results_MinorAP){
  results_MajorAP["wgdStatus"] <- rep(NA, each = length(results_MajorAP$genome_level))
  ploidy <- c(1,2,4,8)
  wgdStatus <- c(0,1,2,3)
  for(i in 1:4){
    results_MajorAP[rowSums(results_MajorAP[,1:39] >= ploidy[i]) >= 20, "wgdStatus"] <- wgdStatus[i]
  }
  loopList <- which(results_MajorAP$wgdStatus == 0)
  for(i in loopList){
    if(rowSums(results_MinorAP[i,1:39] == 0) >= 25 && rowSums(results_MajorAP[i,1:39] > 1) == 0){
      results_MajorAP[i,"wgdStatus"] <- -1
    }
  }
  return(results_MajorAP)
}


#' Get gains, losses and aneuploidyScore (gains + losses) for each sample
#'
#' @description Determines the number of arms which had a gain or a loss and adds data to data frame aneuploidyScoreDF
#'
#' @param results_CopyNumber Data frame with copy number data for each autosome arm of each sample and annotation information
#' @param valNorm Can be either "Genome" or "WGD" and indicates how to normalize the data
#'
#' @return aneuploidyScoreDF
#' @export
#'
#' @examples
#' aneuploidyScores(results_CopyNumber, "Genome")
#' aneuploidyScores(results_CopyNumber, "WGD")
aneuploidyScores <- function(results_CopyNumber, valNorm){
  if(valNorm == "Genome"){
    data_frame <- normalization(results_CopyNumber, TRUE)
  }else if(valNorm == "WGD"){
    data_frame <- wgdNormalization(results_CopyNumber, TRUE)
  }

  y_vector <- rownames(data_frame)

  gains  <- abs(unlist(getCount(data_frame[,1:39], "gain")))
  losses <- abs(unlist(getCount(data_frame[,1:39], "loss")))
  aneuploidyScore <- gains + losses
  gainsLosses <- as.data.frame(cbind(cbind(gains,losses), aneuploidyScore))

  aneuploidyScoreDF <- cbind(cbind(gainsLosses, data_frame$wgdStatus), data_frame$primaryTumorLocation)
  rownames(aneuploidyScoreDF) <- y_vector
  colnames(aneuploidyScoreDF) <- c("gains","losses", "aneuploidyScore", "wgdStatus", "primaryTumorLocation")

  return(aneuploidyScoreDF)
}



#' Get percentage of 0 values in column
#'
#' @description Calculates the percentage of 0 in each column of given data frame and returns
#' list with all percentages of each column in given data frame
#'
#' @param data_frame Can be any given numerical data frame
#'
#' @return tmp_list List with percentages of each column in data frame
#' @export
#'
#' @examples
#' getPercentage(results_CopyNumber)
#' getPercentage(results_MajorAP)
#' getPercentage(results_MinorAP)
getPercentage <- function(data_frame){
  tmp_list <- lapply(c(1:39), function(x){
    totalZero   <- sum(data_frame[,x] == 0)
    percentage  <- (totalZero / length(data_frame[,x])) * 100
    return(percentage)
  })
  return(tmp_list)
}




#' Add primary tumour location
#'
#' @description Merges data frames by row names. Particularly annoCodes and results_Major/Minor or CopyNumber data frame
#'
#' @param data_frame Data frame which has sampleID's as rownames
#' @param annoCodes  Data frame which contains sampleID's and their primary tumour location
#'
#' @return tmpFrame
#' @export
#'
#' @examples
#' addTypeOfCancer(results_CopyNumber, annoCodes)
#' addTypeOfCancer(results_MajorAP,    annoCodes)
#' addTypeOfCancer(results_MinorAP,    annoCodes)
addTypeOfCancer <- function(data_frame, annoCodes){
  tmpFrame <- t(data_frame)
  tmpFrame <- merge(tmpFrame, annoCodes, by= "row.names")
  rownames(tmpFrame) <- tmpFrame[,1]
  tmpFrame <- tmpFrame[,-1]
  return(tmpFrame)
}



#' Get loss or gain count
#'
#' @description Counts the total losses or gains for the data frame depending on given type (loss or gain)
#'
#' @param data_frame Data frame with numercial data
#' @param type Defines if gains or losses need to be counted.
#'
#' @return return_list List with the total losses or gains per data frame column
#' @export
#'
#' @examples
#' getCount(results_CopyNumber, "gain")
#' getCount(results_CopyNumber, "loss")
#' getCount(results_MajorAP,    "gain")
#' getCount(results_MajorAP,    "loss")
#' getCount(results_MinorAP,    "gain")
#' getCount(results_MinorAP,    "loss")
getCount <- function(data_frame, type){
  return_list <- lapply(c(1:length(data_frame$`1p`)), function(x){
    if(type == "gain"){
      if( sum(data_frame[x,] >= 1)){
        total <- sum(data_frame[x,which(data_frame[x,] >=  1)])
      }else{
        total <- 0
      }
    }else if(type == "loss"){
      if(sum(data_frame[x,] <= -1)){
        total <- sum(data_frame[x,which(data_frame[x,] <= -1)])
      }else{
        total <- 0
      }
      return(total)
    }
  })
  return(return_list)
}



#' Data check and repair
#'
#' @description Checks data for NA's in genome_level column and replaces them with most common value of autosomes
#' Identifies positions in df with wrong primaryTumorLocation and changes them to "unknown"
#'
#' @param data_frame Can be either one of the results_ data frames
#'
#' @return data_frame Inserted data frame but without the NA values
#' @export
#'
#' @examples
#' checkDataAndRepair(results_CopyNumber)
#' checkDataAndRepair(results_MajorAP)
#' checkDataAndRepair(results_MinorAP)
checkDataAndRepair <- function(data_frame){
  listNA <- which(is.na(data_frame$genome_level))
  for(x in listNA){
    tmpDF <- table(t(data_frame[x,1:39]))
    data_frame[x, "genome_level"] <- names(tmpDF[which.max(tmpDF)])
  }
  positions <- c(which(grepl("\\d", data_frame$primaryTumorLocation)),
                 which(grepl(".*(male)$", data_frame$primaryTumorLocation)))
  for(i in positions){
    data_frame[i, "primaryTumorLocation"] <- "unknown"
  }
  return(data_frame)
}


#' Get mean of columns
#'
#' @description Calculates mean of every row in the data frame that is given as parameter
#'
#' @param data_frame Can be either one of the results_ data frames
#'
#' @return mean_list List with mean values of the columns of data frame
#' @export
#'
#' @examples
#' getMeans(results_CopyNumber)
#' getMeans(results_MajorAP)
#' getMeans(results_MinorAP)
getMeans <- function(data_frame){
  mean_list <- lapply(c(1:39), function(x){
    chrMean <- mean(data_frame[,x])
  })
  return(mean_list)
}


#' Character string with chromosome arms and genome_level
#'
#' @description Creates character vector with autosome arms
#'
#' @return chromArms A character vector with autosome arms
#' @export
#'
#' @examples
#' getChromArms()
getChromArms <- function(){
  chromArms = c()
  for(i in c(1:22)){
    if(i %!in% c(13,14,15,21,22)){
      chromArms <- c(chromArms, c(paste0(i,"p"), paste0(i,"q")))
    }else{
      chromArms <- c(chromArms, c(paste0(i,"q")))
    }
  }
  return(chromArms)
}



#' Ploidy normalization (genome)
#'
#' @description Normalizes data for the genome level so levels of aneuploidy can be detected in the heatmap plot
#'
#' @param data_frame results_CopyNumber
#' @param valNorm  Boolean which indicates if high values need to be normalized (TRUE = yes, FALSE = no)
#'
#' @return data_frame
#' @export
#'
#' @examples
#' normalization(results_CopyNumber, TRUE)
#' normalization(results_CopyNumber, FALSE)
normalization  <- function(data_frame, valNorm){
  genome_level <- as.numeric(data_frame[,40])
  primaryTumorLocation <- data_frame[,41]
  wgdStatus <- data_frame[,42]

  data_frame <-  data_frame[,1:39]
  data_frame <- (data_frame - genome_level)

  if(valNorm == TRUE){
    data_frame[data_frame >=  1]  <-  1
    data_frame[data_frame <= -1]  <- -1
  }
  data_frame$genome_level <- genome_level
  data_frame$primaryTumorLocation <- primaryTumorLocation
  data_frame$wgdStatus <- wgdStatus

  return(data_frame)
}



#' Ploidy normalization (wgd)
#'
#' @description Normalizes data for the whole genome duplication status so levels of aneuploidy can be detected in the heatmap plot
#'
#' @param data_frame results_CopyNumber data frame
#' @param valNorm Boolean which indicates if high values need to be normalized (TRUE = yes, FALSE = no)
#'
#' @return data_frame
#' @export
#'
#' @examples
#' wgdNormalization(results_CopyNumber, TRUE)
#' wgdNormalization(results_CopyNumber, FALSE)
wgdNormalization  <- function(data_frame, valNorm){
  genome_level <- as.numeric(data_frame[,40])
  primaryTumorLocation <- data_frame[,41]
  wgdStatus  <- data_frame[,42]
  tmpFrame   <- data_frame
  data_frame <- data_frame[,1:39]

  data_frame[which(tmpFrame[,42] ==-1),] <- data_frame[which(tmpFrame[,42] ==-1),]  - 1
  data_frame[which(tmpFrame[,42] == 0),] <- data_frame[which(tmpFrame[,42] == 0),]  - 2
  data_frame[which(tmpFrame[,42] == 1),] <- data_frame[which(tmpFrame[,42] == 1),]  - 4
  data_frame[which(tmpFrame[,42] == 2),] <- data_frame[which(tmpFrame[,42] == 2),]  - 8

  if(valNorm == TRUE){
    data_frame[data_frame >=  1]  <-  1
    data_frame[data_frame <= -1]  <- -1
  }
  data_frame$genome_level <- genome_level
  data_frame$wgdStatus    <- wgdStatus
  data_frame$primaryTumorLocation <- primaryTumorLocation
  return(data_frame)
}



#' Clean data
#'
#' @description Filters out all cancertypes with a prevalence less or equal than 10 and the unknown or wrongly classified
#'
#' @param data_frame Can be either one of results_ data frames
#'
#' @return data_frame
#' @export
#'
#' @examples
#' cleanData(results_CopyNumber)
#' cleanData(results_MajorAP)
#' cleanData(results_MinorAP)
cleanData <- function(data_frame){
  data_frame$primaryTumorLocation  <- sapply(data_frame$primaryTumorLocation, tolower)
  uniqueTypes <- unique(data_frame$primaryTumorLocation)

  total_list  <- lapply(uniqueTypes, function(type){
    return(sum(data_frame$primaryTumorLocation == type))
  })
  prevalence_df <- as.data.frame(cbind(uniqueTypes, total_list))
  rmCancerTypes <- c(unlist(prevalence_df$uniqueTypes[prevalence_df$total_list <= 10]), "unknown")

  lapply(rmCancerTypes, function(tp){
    data_frame <<- data_frame[!(data_frame$primaryTumorLocation == tp),]
  })
  return(data_frame)
}



#' Get mutational load annotations
#'
#' @description Constructs one mutational load data frame for annotations per sample of three single mutational load data frames
#'
#' @param mut_load_df Data frame which contains mutational load data per sample
#' @param clonal_mut_load_df Data frame which contains clonal and subclonal mutational load data per sample
#' @param sign_contribution_df Data frame which contains the signature contribution for each signature per sample
#'
#' @return mut_load_anno Data frame with mutational load annotation per sample
#' @export
#'
#' @examples
#' getMutloadAnnotations(mut_load_df, clonal_mut_load_df, sign_contribution_df)
getMutloadAnnotations <- function(mut_load_df, clonal_mut_load_df, sign_contribution_df){
  mut_load_anno <- cbind(clonal_mut_load_df, mut_load_df$mut_load)
  colnames(mut_load_anno) <- c(colnames(clonal_mut_load_df), "mut_load_total")
  sign_anno <- sign_contribution_df[,31:60]
  mut_load_anno <- cbind(mut_load_anno, sign_anno)
  colnames(mut_load_anno) <- c(colnames(clonal_mut_load_df), "mut_load_total", colnames(sign_anno))
  return(mut_load_anno)
}


#' Add mutational load data to aneuploidyScoreDF
#'
#' @description Adds mutational load data to aneuploidyScoreDF
#'
#' @param aneuploidyScoreDF Data frame which contains gains, losses and aneuploidy score (gains + losses) per sample
#' @param mut_load_df Data frame which contains mutational load data per sample
#' @param clonal_mut_load_df Data frame which contains clonal and subclonal mutational load data per sample
#'
#' @return aneuploidyScoreDF
#' @export
#'
#' @examples
#' getAneuploidyDF(aneuploidyScoreDF, mut_load_df, clonal_mut_load_df)
getAneuploidyDF <- function(aneuploidyScoreDF, mut_load_df , clonal_mut_load_df){
  aneuploidyScoreDF <- cbind(aneuploidyScoreDF, mut_load_df[,"mut_load"][match(rownames(aneuploidyScoreDF), rownames(mut_load_df))])
  aneuploidyScoreDF <- cbind(aneuploidyScoreDF, clonal_mut_load_df[,1][match(rownames(aneuploidyScoreDF)  , rownames(clonal_mut_load_df))])
  aneuploidyScoreDF <- cbind(aneuploidyScoreDF, clonal_mut_load_df[,2][match(rownames(aneuploidyScoreDF)  , rownames(clonal_mut_load_df))])

  colnames(aneuploidyScoreDF) <- c("gains", "losses", "aneuploidyScore", "wgdStatus","primaryTumorLocation",
                                   "mut_load", "clonal_mut_load", "subclonal_mut_load")
  return(aneuploidyScoreDF)
}

