% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_alleles.R
\name{check_alleles}
\alias{check_alleles}
\title{check alleles}
\usage{
check_alleles(ref, alt, r.ref, r.alt)
}
\arguments{
\item{ref}{ref allele for target data}

\item{alt}{allele for target data}

\item{r.ref}{ref allele for reference data}

\item{r.alt}{alt allele for reference data}
}
\value{
genotpe status based on allele labels
}
\description{
Checks if alleles match, are strand flipped, switched (ref/alt swap)
}
\examples{
check_alleles("A", "T", "T", "A")
}
