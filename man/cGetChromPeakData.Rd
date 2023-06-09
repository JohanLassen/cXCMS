% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parallel_peak_filling.R
\name{cGetChromPeakData}
\alias{cGetChromPeakData}
\title{Intensive peak extraction step performed file-wise (parallel in cluster workflows)}
\usage{
cGetChromPeakData(prepared_data, index, output_dir = "./tmp/integration_tmp/")
}
\arguments{
\item{prepared_data}{as returned by cXCMS::cPrepFillChromPeaks()}

\item{index}{the file index in sample_info}

\item{output_dir}{the tmp directory for storage of extracted peaks}
}
\value{
a res object for the final peak integration.
}
\description{
Intensive peak extraction step performed file-wise (parallel in cluster workflows)
}
