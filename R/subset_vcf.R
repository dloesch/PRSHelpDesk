#'Subset vcf file at given position
#'
#' Subsets vcf file via system call for bcftools. File needs to be indexed
#' and bcftools needs to be installed and in path
#' @importFrom utils read.table
#' @param input.vcf filepath to vcf file
#' @param chr chromosome interger
#' @param pos physical poisition
#' @param snp snp ID. Defaults to NA
#' @param build genome build. Defaults to hg19
#' @return data.frame of vcf at given posiiton
#' @export
#' @examples
#' \dontrun{
#' subset_vcf("data.vcf.gz", chr=2, pos=756111)
#' }

subset_vcf <- function(input.vcf, chr, pos, snp=NA, build="hg19"){
  if(tolower(build) != "hg19") chr <- paste0("chr", chr)
  p <- pipe(paste("bcftools view -r", paste0(chr, ":", pos), input.vcf), "rt")
  v <- read.table(p, colClasses = c("character"))
  close(p)
  if(!is.na(snp)){
    v <- v[v[,3] == snp,]
  }

  return(v)
}
