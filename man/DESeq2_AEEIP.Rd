% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DESeq2_AEEIP.R
\name{DESeq2_AEEIP}
\alias{DESeq2_AEEIP}
\title{DESeq2 with correct bias}
\usage{
DESeq2_AEEIP(AEEIP_test_final_list, input, IP, times = 100, seed = 1)
}
\arguments{
\item{AEEIP_test_final_list}{list of AEEIP_test_final result of different replicates}

\item{input}{input read counts of different replicates}

\item{IP}{IP read counts of different replicates}

\item{times}{sample times}

\item{seed}{sample seed}
}
\value{
DESeq2 results
}
\description{
DESeq2 with correct bias
}
