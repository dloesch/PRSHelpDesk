% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annoLAI.R
\name{annotLAI}
\alias{annotLAI}
\title{Annotate VCF with LAI}
\usage{
annotLAI(
  vcf.file,
  msp.file,
  chr,
  out.prefix = "output",
  compress = TRUE,
  bgzip = "bgzip"
)
}
\arguments{
\item{vcf.file}{path to vcf file. Note: vcf file should just contain phased genotypes}

\item{msp.file}{path to msp file from rfmix}

\item{chr}{chromsoome number}

\item{out.prefix}{prefix for output vcf}

\item{compress}{compress file? Defaults to true}

\item{bgzip}{path to utility to compress file. Defaults to "bgzip"}
}
\description{
Annotates VCF with local ancestry information by appending results from msp file (LAI) with VCF genotpe records.
Can be paired with laiPRS() to estimate a local-ancestry aware PRS
}
