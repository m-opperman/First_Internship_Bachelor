###########################################################################################################################################
# Author:   M. Opperman                                                                                                                   #
# Date:     09-10-2018                                                                                                                    #
# Version:  3.0                                                                                                                           #
#                                                                                                                                         #
# Function: Contains core functions for reading files and determination of autosome ploidy levels.                                        #                                                                                                                                         #
###########################################################################################################################################



#' Reading .cnv file
#'
#' @description Reads file and makes data frame of data (sample_data)
#'
#' @param filename Character string of the file that need to be read.
#'
#' @return sample_data
#' @export
#'
#' @examples
#' readFile("example.cnv")
readFile <- function(filename){
  #  filename <- list.files(pattern = filename, path = "/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data", recursive = TRUE, full.names = T)
  first_line  <- readLines(filename, n=1)
  sample_data <- read.table(filename, header=FALSE, sep= "\t")
  names(sample_data) <- as.vector(unlist(strsplit(first_line[length(first_line)],"\t")))
  return(sample_data)
}



#' Add Minor and Major allele ploidy levels
#'
#' @description Adds major and minor allele ploidy data to the data frame
#'
#' @param sample_data Data frame with all information from the .cnv file
#'
#' @return sample_data
#' @export
#'
#' @examples
#' addMajorMinorAllele(sample_data)
addMajorMinorAllele <- function(sample_data){
  sample_data$copyNumber[sample_data$copyNumber < 0] <- 0
  sample_data["unboundMinorAllele"] <- (1 - sample_data$actualBAF) * sample_data$copyNumber
  sample_data$unboundMinorAllele[sample_data$unboundMinorAllele < 0.25 ] <- 0
  sample_data["majorAllele"] <- sample_data$copyNumber - sample_data$unboundMinorAllele
  sample_data$majorAllele[sample_data$majorAllele < 0 ] <- 0
  return(sample_data)
}


#' Determine ploidy levels
#'
#' @description Calculates per chromosome arm the average major and minor ploidy (weighted)
#'
#' @param sample_data Data frame with all information from the .cnv file
#' @param filename Character string of the file that need to be read.
#'
#' @return list with major and minor ploidy levels per autosome arm
#' @export
#'
#' @examples
#' chromArmPloidy(sample_data, "example.cnv")
chromArmPloidy <- function(sample_data, filename){
  majorPloidy_list <- c()
  minorPloidy_list <- c()
  lapply(c(1:22), function(c){
    if(c %!in% c(13,14,15,21,22)){
      subset_data <- subset(sample_data, sample_data$`#chromosome` == c)
      loc_centromeer <- which(subset_data$segmentEndSupport == "CENTROMERE")


      #Calculates everything for P arm
      subset_p  <- subset_data[1:loc_centromeer,]
      df_subset <- majorGroup(subset_p)

      ploidyFrameP <- createPloidyFrame(df_subset, c, "p")
      return_listP <- aneuPloidy(ploidyFrameP, c, "p", filename)

      majorPloidy_list[length(majorPloidy_list)+1] <<- return_listP[1]
      minorPloidy_list[length(minorPloidy_list)+1] <<- return_listP[2]


      #Calculates everything for Q arm
      subset_q   <- subset_data[(loc_centromeer+1):length(subset_data$start),]
      df_subset2 <- majorGroup(subset_q)

      ploidyFrameQ <- createPloidyFrame(df_subset2,c, "q")
      return_listQ <- aneuPloidy(ploidyFrameQ, c, "q", filename)

      majorPloidy_list[length(majorPloidy_list)+1] <<- return_listQ[1]
      minorPloidy_list[length(minorPloidy_list)+1] <<- return_listQ[2]

    }else{

      #Calculates everything for the autosomes that are treated as one armed
      subset_data <- subset(sample_data, sample_data$`#chromosome` == c)
      df_subset   <- majorGroup(subset_data)

      ploidyFrame <- createPloidyFrame(df_subset, c, "q")
      return_list <- aneuPloidy(ploidyFrame, c, "q", filename)

      majorPloidy_list[length(majorPloidy_list)+1] <<- return_list[1]
      minorPloidy_list[length(minorPloidy_list)+1] <<- return_list[2]
    }
  })
  return(list(majorPloidy_list, minorPloidy_list))
}


#' Optimalize the major and minor ploidy levels
#'
#' @description Determines overall MajorAP level of the genome and changes the NA values accordingly with previous stores doubt_df data.
#'
#' @param ploidy_list list of ploidy levels
#' @param filename Character string of the file that need to be read.
#' @param tp Character string that defines if major or minor doubt cases get normalized
#'
#' @return ploidy_list
#' @export
#'
#' @examples
#' optimalization(ploidy_list, "example.cnv", "major")
optimalization <- function(ploidy_list, filename, tp){
  ploidy_list  <- as.numeric(ploidy_list)
  common_value <- names(sort(summary(as.factor(ploidy_list)), decreasing=TRUE, na.last = TRUE))[1:2]
  genome_level <- common_value[1]
  ploidy_list  <- append(ploidy_list, values = genome_level)
  chromArms    <- getChromArms()
  tmp_df <- as.data.frame(cbind(chromArms, ploidy_list), stringsAsFactors = FALSE)

  doubt_df <- doubt_df[(doubt_df$file == filename),]
  doubt_df <- doubt_df[(doubt_df$tp == tp),]
  if(length(doubt_df$chromArm) == 0){
    return(ploidy_list)
  }else{
    for(i in 1:nrow(doubt_df)){
      row <- as.list(doubt_df[i,])
      if(row$highest == genome_level){
        tmp_df[which(tmp_df$chromArms == row$chromArm), "ploidy_list"] <- genome_level
      }else if(row$second == genome_level){
        tmp_df[which(tmp_df$chromArms == row$chromArm), "ploidy_list"] <- genome_level
      }else{
        tmp_df[which(tmp_df$chromArms == row$chromArm), "ploidy_list"] <- row$highest
      }
      tmp_df$ploidy_list <- as.numeric(tmp_df$ploidy_list)
    }
    ploidy_list <- tmp_df$ploidy_list
  }
  return(ploidy_list)
}
