% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bar_boxPlots.R
\name{violinAneuploidy}
\alias{violinAneuploidy}
\title{Distribution of gains, losses and aneuploidyScore per WGD status (Violin plot)}
\usage{
violinAneuploidy(aneuploidyScoreDF_wgd)
}
\arguments{
\item{aneuploidyScoreDF_wgd}{Data frame which contains gains, losses and primary wgdStatus for each sample.}
}
\value{
None
}
\description{
Plots distribution of gains and losses per wgd Status and primaryTumorLocation
}
\examples{
violinAneuploidy(aneuploidyScoreDF_wgd)
violinAneuploidy(aneuploidyScoreDF)
}
