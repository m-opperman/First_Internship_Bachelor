###########################################################################################################################################
# Author:   M. Opperman                                                                                                                   #
# Date:     09-10-2018                                                                                                                    #
# Version:  3.0                                                                                                                           #
#                                                                                                                                         #
# Function: Estimates the ploidy level per autosome arm for each sample.                                                                  #
#           Input = .cnv file (with derived data from WGS of cancer patient) with for each segment the following info:                    #
#           #chromosome	,start	,end	,copyNumber	,bafCount	,observedBAF	,actualBAF	,segmentStartSupport	,segmentEndSupport	,method   #
#           This information is used to calculate the minor and major allele ploidy of each segment and finally for each autosome arm.    #
#                                                                                                                                         #
# Output:   Gives three .tsv files as output with per sample the major, minor or copy number of each autosome arm:                        #
#           Results_MajorAP.tsv  /  Results_MinorAP.tsv  /  Results_CopyNumber.tsv                                                        #
###########################################################################################################################################

#Necessary scripts
source("coreFunctions.R")
source("supportFunctions.R")

#Imports
library("data.table")

#Global settings, functions and variables
options(stringsAsFactors = FALSE)
options(scipen = 999)


'%!in%' <- function(x,y)!('%in%'(x,y))
doubt_df <- data.frame(file = character(),chromArm = character(), highest= integer(), second= integer(), tp = character(),
                       stringsAsFactors = FALSE)



#' Main
#'
#' @description Main function which calls all necessary functions for determination of ploidy levels per autosome
#'
#' @return nothing but does write PDF files with all plots of a file and a summary pdf of all files.
#' @export
#'
#' @examples
ploidyEstimator <- function(path){
  files <- list.files(path = path,  pattern = "\\.purple.cnv$", recursive = TRUE, full.names = TRUE)

  chromArms <- getChromArms()
  df_majorPloidy <- data.frame(chromArms)
  df_minorPloidy <- data.frame(chromArms)
  df_copyNumber  <- data.frame(chromArms)

  return_list <- lapply(files, function(file){

    sample_data <- readFile(file)
    sample_data <- addMajorMinorAllele(sample_data)
    return_list <- chromArmPloidy(sample_data, file)

    majorPloidyList <- optimalization(unlist(return_list[1]), file, "major")
    minorPloidyList <- optimalization(unlist(return_list[2]), file, "minor")

    df_majorPloidy[file] <<- as.numeric(majorPloidyList)
    df_minorPloidy[file] <<- as.numeric(minorPloidyList)
    df_copyNumber[file]  <<- as.numeric(as.character(majorPloidyList)) + as.numeric(as.character(minorPloidyList))

  })

  write.table(df_majorPloidy, file = "Results_MajorAP.tsv"    , sep="\t", col.names = TRUE, row.names = FALSE)
  write.table(df_minorPloidy, file = "Results_MinorAP.tsv"    , sep="\t", col.names = TRUE, row.names = FALSE)
  write.table(df_copyNumber , file = "Results_CopyNumber.tsv" , sep="\t", col.names = TRUE, row.names = FALSE)
}


ploidyEstimator(path = "C:/Users/Username/Desktop/R/Packages/ploidyEstimator/cnv_files/")
