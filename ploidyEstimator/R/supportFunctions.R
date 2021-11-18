###########################################################################################################################################
# Author:   M. Opperman                                                                                                                   #
# Date:     09-10-2018                                                                                                                    #
# Version:  3.0                                                                                                                           #
#                                                                                                                                         #
# Function: Contains support functions for determining the ploidy levels of each autosome arm.                                            #                                                                                                                                         #
###########################################################################################################################################


#' Determine groups
#'
#' @description Assigns each segment to the right group
#'
#' @param data_frame Data frame with all information from the .cnv file
#'
#' @return data_frame
#' @export
#'
#' @examples
#' majorGroup(data_frame)
majorGroup <- function(data_frame){
  data_frame["majorGroup"] <- cut(data_frame$majorAllele,
                                 breaks= c(-1.5,-0.5,0.5,c(1*(1.5:20))), right = FALSE,
                                 include.lowest = TRUE, labels = c(-1:19))

  data_frame["group"] <- rleid(data_frame$majorGroup)
  data_frame$majorAllele[data_frame$majorAllele <= 0 ] <- NA
  data_frame["mapWeight"] <- rep(1, length(data_frame[,1]))
  data_frame$mapWeight[is.na(data_frame$mapWeight)] <- 1

  return(data_frame)
}


#' Combine segments based on ploidy level
#'
#' @description Creates new table with combined segments if it's possible
#'
#' @param arm_subset Data frame with segment information of one chromosome arm
#' @param chromosome Integer that indicates which chromosome is represented in arm_subset
#' @param arm Character that indicates which chromosome arm is represented in arm_subset
#'
#' @return ploidyFrame
#' @export
#'
#' @examples
#' smoothing(arm_subset, 13 , "q")
createPloidyFrame <- function(arm_subset, chromosome, arm){
  grouplist <- unique(arm_subset$group)

  tmp_list  <- lapply(1:length(unique(arm_subset$group)), function(i){
    group_set <- subset(arm_subset, arm_subset$group == grouplist[i])
    start <- group_set$start[1]
    end   <- group_set$end[length(group_set$end)]

    segment_size <- (group_set$end - group_set$start) + 1
    weight <-  segment_size * group_set$mapWeight

    averageMAP <- group_set$majorAllele * weight
    averageMAP <- sum(averageMAP) / sum(weight)

    averageMinor <- group_set$unboundMinorAllele * weight
    averageMinor <- sum(averageMinor) / sum(weight)

    averageCN <- group_set$copyNumber * weight
    averageCN <- sum(averageCN) / sum(weight)

    if(!is.na(averageMinor) && averageMinor <= 0.25){
      loh <- TRUE
    }else{
      loh <- FALSE
    }
    tmp_vector <- c(-1:19)
    return(c(paste0(chromosome,arm), as.integer(start), as.integer(end), as.numeric(averageMAP), as.numeric(averageMinor),
             as.numeric(averageCN), as.integer(tmp_vector[group_set$majorGroup[1]]), as.logical(loh)))
  })

  ploidyFrame <- as.data.frame(do.call(rbind, tmp_list), stringsAsFactors = FALSE)
  colnames(ploidyFrame) <- c("chromArm", "start", "end", "averageMajorAP", "averageMinorAP", "averageCN", "majorAPcategory", "lohStatus")
  ploidyFrame["minorAPcategory"] <- round(as.numeric(ploidyFrame$averageMinorAP))

  return(ploidyFrame)
}



#' Determine ploidy level of chromosome arm
#'
#' @description Determines ploidy level of the autosome arms
#'
#' @param ploidyFrame Data frame with segment information of one chromosome arm
#' @param chromosome Integer that indicates which chromosome is represented in arm_subset
#' @param arm Character that indicates which chromosome arm is represented in arm_subset
#' @param filename Character string of the file that need to be read.
#'
#' @return list with major and minor allele ploidy levels
#' @export
#'
#' @examples
#' aneuPloidy(ploidyFrame, 12 , "p", "example.cnv")
aneuPloidy <- function(ploidyFrame, chromosome, arm, filename){
  ploidyFrame["segmentSize"] <- (as.numeric(ploidyFrame$end) - as.numeric(ploidyFrame$start))

  major_list <-  aggregate(segmentSize~majorAPcategory, data = ploidyFrame, FUN = sum)
  major_list["percentage"]   <- (major_list$segmentSize / sum(major_list$segmentSize)) * 100

  minor_list <-  aggregate(segmentSize~minorAPcategory, data = ploidyFrame, FUN = sum)
  minor_list["percentage"]   <- (minor_list$segmentSize / sum(minor_list$segmentSize)) * 100

  tmpList  <- list(major_list, minor_list)
  category <- c("major", "minor")

  counter <- 1
  ploidyList <- lapply(tmpList, function(subList){

    column <- paste0(category[counter], "APcategory")
    if(length(order(subList$percentage)) > 1){
      if(max(subList$percentage) >= 50  && (subList[which.max(subList$percentage), "percentage"] - sort(subList$percentage, decreasing = TRUE)[2]) > 10){
        ploidyLevel <- subList[which.max(subList$percentage), column]
      }else{
        ploidyLevel <- NA
        doubt_df[nrow(doubt_df)+1,] <<- c(filename, paste0(chromosome, arm),  subList[which.max(subList$percentage), column],
                                          subList[which(subList$percentage == (sort(subList$percentage,decreasing=TRUE)[2])), column], category[counter])
      }
    }else{
      ploidyLevel <- subList[which.max(subList$percentage), column]
    }
    counter <<- counter + 1
    return(ploidyLevel)
  })

  return(list(as.numeric(unlist(ploidyList[1])), as.numeric(as.character(unlist(ploidyList[2])))))
}


#' Get autosome arm string
#'
#' @description Creates vector containing autosome arms
#'
#' @return chromArms
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
  chromArms <- append(chromArms, "genome_level")
  return(chromArms)
}
