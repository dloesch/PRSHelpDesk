#' rename SNPs
#'
#' Takes advantage of the SNP matching routine in constructing mgData class object to write out SNP ids from target and gwas stats
#' for updating target SNP ids. Can also write out a filtered vcf file with new ids, if desired.
#'
#' @importFrom data.table fwrite fread
#' @param mgData mgData object
#' @param CGAT keep CG/AT sites. defaults to FALSE
#' @param vcf.file file path to vcf file. Defaults to NA. If provided, will read in and write out vcf with updated ids.
#' @param prefix output prefix. Defaults to "output"
#' @param compress Compress VCF? defaults to true
#' @param bgzip Path to bgzip or gzip
#' @export
#'
renameSNPs <- function(mgData, CGAT=FALSE, vcf.file=NA, prefix="output", compress=TRUE, bgzip="bgzip"){

  x <- mgData

  if(CGAT == FALSE){
    x <- subset(x, x$cgat == FALSE)
  }

  #filter
  x <- subset(x, x$filter == "PASS")

  snps <- data.frame(x$id1, x$id2, stringsAsFactors = FALSE)

  if(!is.na(vcf.file)){
    vcf <- read_vcf(vcf.file)
    vcf <- vcf[vcf[,3] %in% x$id1,]
    vcf[,3] <- x$id2
    write_vcf(vcf, paste0(prefix, ".vcf"))

    if(compress == "TRUE"){
      system2(bgzip, paste("-f", paste0(prefix, ".vcf")))
    }

  }

  return(snps)

}
