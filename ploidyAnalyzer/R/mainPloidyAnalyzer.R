###########################################################################################################################################
# Author:   M. Opperman                                                                                                                   #
# Date:     20-11-2018                                                                                                                    #
# Version:  2.0                                                                                                                           #
#                                                                                                                                         #
# Function: Main function for analyses of the ploidy levels.                                                                              #
#           Input = three .tsv files      (Results_MajorAP, Results_MinorAP, Results_CopyNumber)                                          #
#                   four  .txt/.csv files (mut_matrix_and_mut_load, mut_matrix_andmut_load_clonal_subclonal,                              #
#                                          mut_matrix_and_sign_contribution, HMF_Annotations)                                             #
#           The "Results_" .tsv files contain the calculated ploidy levels per autosomearm of each cancer sample                          #
#           The annotation files contain information about: primary tumor location, mutational load and signature contribution.           #                                             
#                                                                                                                                         #
# Output:   Gives one pdf file with all the diagnostic plots in it.                                                                       #
###########################################################################################################################################

# #R script which are needed for the analysis
# source("readFunctions.R")
# source("supportFunctions.R")
# source("heatmapPlots.R")
# source("bar_boxPlots.R")
# source("pre_proces_heatmaps.R")
#
# #Necessary libraries
# library("lattice")
# library("stringr")
# library("ggplot2")
# library("heatmap.plus")
# library("gplots")
# library("reshape2")

#Global settings and functions
options(stringsAsFactors = FALSE)
'%!in%' <- function(x,y)!('%in%'(x,y))


#' Clustering function
#'
#' @description Average cluster function for clustering the unsupervised heatmaps
#'
#' @param data_frame The data frame which needs to be clustered
#'
#' @return None
#' @export
#'
#' @examples
#' hclust()
hclust2 <- function(data_frame){
  hclust(data_frame, method = "average")
}

#Input and output path for writing and reading files from other directories
outputPath <- "~/Documents/Packages/ploidyAnalyzer/output/"
#outputPath <- "C:/Users/melan/Desktop/Stage/R/Packages/PloidyAnalyzer/output/"
inputPath  <- "~/Documents/Packages/ploidyAnalyzer/input/"
#inputPath  <- "C:/Users/melan/Desktop/Stage/R/Packages/PloidyAnalyzer/input/"


#' Main function of ploidyAnalyzer
#'
#' @description Calls every function in the right order and produces .pdf and .png files with all the plots
#'
#' @param results_MajorAP        .tsv file which contains major allele ploidy levels for each autosome of each sample and their genome_level
#' @param results_MinorAP        .tsv file which contains minor allele ploidy levels for each autosome of each sample and their genome_level
#' @param results_CopyNumber     .tsv file which contains copy number levels for each autosome of each sample and their genome_level
#' @param mutational_load        .csv file which contains mutational load data for each sample
#' @param clonal_mutational_load .csv file which contains clonal and subclonal mutational load data for each sample
#' @param signature_contribution .csv file which contains the (relative) signature contribution for each sample
#' @param msi_patients           .txt file which contains the sampleID of patients which are marked as MSI
#' @param annotations            .txt file which contains the sampleID's and their primary tumour location
#'
#' @return None
#' @export
#'
#' @examples
#' ploidyAnalyzer("Results_MajorAP.tsv","Results_MinorAP.tsv","Results_CopyNumber.tsv","mut_matrix_and_mut_load.csv","mut_matrix_and_mut_load_clonal_subclonal.csv","mut_matrix_and_sign_contribution.csv","msi_patients.txt","HMF_Annotations.txt")
ploidyAnalyzer <- function(results_MajorAP, results_MinorAP, results_CopyNumber, mutational_load, clonal_mutational_load, signature_contribution, msi_patients, annotations){
  #Reads all annotation files
  mut_load_df           <- as.data.frame(read.table(file = paste0(inputPath, mutational_load)       ,sep = "\t", header = TRUE, row.names = 1))
  clonal_mut_load_df    <- as.data.frame(read.table(file = paste0(inputPath, clonal_mutational_load),sep = "\t", header = TRUE, row.names = 1))
  sign_contribution_df  <- as.data.frame(read.table(file = paste0(inputPath, signature_contribution),sep = "\t", header = TRUE, row.names = 1))
  msiPatients           <- scan(paste0(inputPath, msi_patients), character(), quote = "\n")
  annoCodes             <- readAnnoCodes(paste0(inputPath, annotations))

  #Reads all results .tsv files from ploidy Estimator script and adds Whole genome duplication status
  results_MajorAP    <- readFile(paste0(inputPath, results_MajorAP)   , annoCodes)
  results_MinorAP    <- readFile(paste0(inputPath, results_MinorAP)   , annoCodes)
  results_CopyNumber <- readFile(paste0(inputPath, results_CopyNumber), annoCodes)

  results_MajorAP    <- addSampleState(results_MajorAP, results_MinorAP)
  results_MinorAP["wgdStatus"]    <- results_MajorAP$wgdStatus
  results_CopyNumber["wgdStatus"] <- results_MajorAP$wgdStatus

  #Creates aneuploidyScoreDF with all absolute gains and losses per sample
  aneuploidyScoreDF_genome <- aneuploidyScores(results_CopyNumber, "Genome")
  aneuploidyScoreDF_wgd    <- aneuploidyScores(results_CopyNumber, "WGD")

  #Creates heatmaps.pdf and calls all functions for supervised clustered heatmaps
  pdf(paste0(outputPath,"heatmaps.pdf"), width = 20, height = 15)
  #pdf("heatmaps.pdf", width = 20, height = 15)

  plotCopyNumber(results_CopyNumber, "CopyNumber")
  plotLOH(results_MinorAP)

  dev.off()

  aneuploidyScoreDF_wgd <- getAneuploidyDF(aneuploidyScoreDF_wgd  ,  mut_load_df ,  clonal_mut_load_df)
  mut_load_anno         <- getMutloadAnnotations(mut_load_df, clonal_mut_load_df, sign_contribution_df)

  #Creates bar_boxplots.pdf and calls all functions for bar, box and violin plots
  pdf(paste0(outputPath,"bar_boxplots.pdf"), width = 15, height = 8)
  #pdf("bar_boxplots.pdf", width = 15, height = 8)

  duplicated_arms_plot(results_MajorAP)
  cancerTypes(results_CopyNumber)
  outliers(results_CopyNumber, 10)
  lohDistribution(results_MinorAP)
  violinAneuploidy(aneuploidyScoreDF_wgd)

  dev.off()

  #Functions for creating .png files with unsupervised heatmaps
  lohUnsupervisedHeatmap(results_MinorAP, mut_load_anno)
  preProcesUnsupervisedHeatmap(results_CopyNumber, mut_load_anno, "WGD")
  preProcesUnsupervisedHeatmap(results_CopyNumber, mut_load_anno, "NonWGD")
  preProcesUnsupervisedHeatmap(results_CopyNumber, mut_load_anno, "All")

}

ploidyAnalyzer("Results_MajorAP.tsv","Results_MinorAP.tsv","Results_CopyNumber.tsv","mut_matrix_and_mut_load.csv",
               "mut_matrix_and_mut_load_clonal_subclonal.csv","mut_matrix_and_sign_contribution.csv","msi_patients.txt",
               "HMF_Annotations.txt")
#main()
