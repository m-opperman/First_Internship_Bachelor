% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mainPloidyAnalyzer.R
\name{ploidyAnalyzer}
\alias{ploidyAnalyzer}
\title{Main function of ploidyAnalyzer}
\usage{
ploidyAnalyzer(results_MajorAP, results_MinorAP, results_CopyNumber,
  mutational_load, clonal_mutational_load, signature_contribution,
  msi_patients, annotations)
}
\arguments{
\item{results_MajorAP}{.tsv file which contains major allele ploidy levels for each autosome of each sample and their genome_level}

\item{results_MinorAP}{.tsv file which contains minor allele ploidy levels for each autosome of each sample and their genome_level}

\item{results_CopyNumber}{.tsv file which contains copy number levels for each autosome of each sample and their genome_level}

\item{mutational_load}{.csv file which contains mutational load data for each sample}

\item{clonal_mutational_load}{.csv file which contains clonal and subclonal mutational load data for each sample}

\item{signature_contribution}{.csv file which contains the (relative) signature contribution for each sample}

\item{msi_patients}{.txt file which contains the sampleID of patients which are marked as MSI}

\item{annotations}{.txt file which contains the sampleID's and their primary tumour location}
}
\value{
None
}
\description{
Calls every function in the right order and produces .pdf and .png files with all the plots
}
\examples{
ploidyAnalyzer("Results_MajorAP.tsv","Results_MinorAP.tsv","Results_CopyNumber.tsv","mut_matrix_and_mut_load.csv","mut_matrix_and_mut_load_clonal_subclonal.csv","mut_matrix_and_sign_contribution.csv","msi_patients.txt","HMF_Annotations.txt")
}
