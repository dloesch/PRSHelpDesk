% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parseVCF.R
\name{parseVCF}
\alias{parseVCF}
\title{parse genotype files}
\usage{
parseVCF(
  geno.file,
  build = "hg19",
  geno = FALSE,
  both = FALSE,
  af = FALSE,
  plink = "plink",
  bcftools = "bcftools"
)
}
\arguments{
\item{geno.file}{file path to genotype data}

\item{build}{genome build. Either hg19 or hg38. Defaults to hg19}

\item{geno}{Read in genotype data? Defaults to FALSE}

\item{both}{Read in full data? Defaults to FALSE}

\item{af}{Calculate allele frequencies? Deafuls to FALSE}

\item{plink}{path to plink}

\item{bcftools}{path to bcftools}
}
\description{
parses genotype files, converting either to class gData if by-site information is wanted
or class vcfData if genotype data is wanted
}
