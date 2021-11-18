###########################################################################################################################################
# Author:   M. Opperman                                                                                                                   #
# Date:     20-11-2018                                                                                                                    #
# Version:  2.0                                                                                                                           #
#                                                                                                                                         #
# Function: Contains the necessary read functions to analyse the ploidy levels of cancer samples.                                         #
#           These functions are being called from the main_PA.R script                                                                    #
###########################################################################################################################################


#' Read annotation file with sampleID's and primary tumour locations
#'
#' @description Reads file with annotations (primaryTumorLocation of each sample) and turns it into a dataframe with sampleId and primaryTumorLocation
#'
#' @param filename .txt file which contains sampleID column and primary tumor location column.
#'
#' @return annoCodes
#' @export
#'
#' @examples
#' readAnnoCodes(filename)
readAnnoCodes <- function(filename){
  first_line  <- readLines( filename, n=1)
  annoCodes   <- read.table(filename, header = TRUE, sep= " ", row.names = 1)
  names(annoCodes) <- as.vector(unlist(strsplit(first_line[length(first_line)]," ")))[2]
  annoCodes$primaryTumorLocation <- sapply(annoCodes$primaryTumorLocation, tolower)

  return(annoCodes)
}



#' Read .tsv files
#'
#' @description Reads results .tsv files and makes a data frame
#'
#' @param filename .tsv results file from ploidyEstimator script
#' @param annoCodes dataframe with sampleID and primary tumor location data
#'
#' @return data_frame
#' @export
#'
#' @examples
#' readFile(filename, annoCodes)
readFile <- function(filename, annoCodes){
  first_line    <- readLines(filename, n=1)
  column_vector <- str_extract_all(str_extract_all(first_line, "[^/]*(.purple.cnv)", simplify = TRUE), ".*[^(.purple.cnv)]",  simplify = TRUE)
  data_frame    <- read.table(filename, header=TRUE, sep= "\t", row.names = 1)
  names(data_frame) <- column_vector
  data_frame <- addTypeOfCancer(data_frame, annoCodes)
  data_frame <- checkDataAndRepair(data_frame)
  data_frame <- cleanData(data_frame)
  return(data_frame)
}






