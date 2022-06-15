#' Write VCF
#'
#' Wites VCF from vcf in memory
#' @importFrom data.table fwrite
#' @param vcf data.frame of vcf file in memory
#' @param out.file name of new vcf
#' @return new vcf file
#' @export
#' @examples
#' \dontrun{
#' write_vcf2(vcf, "data.vcf")
#' }
write_vcf <- function(vcf, out.file){

  vcf$FILTER <- "PASS"

  #prepare output vcf
  if(length(vcf$INFO[1][grep("AF", vcf$INFO[1])]) == 1){
    output.header <- c("##fileformat=VCFv4.2",
                       paste0("##filedate=", format(Sys.time(), "%Y%m%d")),
                       "##source=PRStools",
                       "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>",
                       "##INFO=<ID=AF,Number=A,Type=Float,Description='Allele frequency'>",
                       paste(colnames(vcf), collapse="\t"))
  }else{
    output.header <- c("##fileformat=VCFv4.2",
                       paste0("##filedate=", format(Sys.time(), "%Y%m%d")),
                       "##source=PRStools",
                       "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>",
                       paste(colnames(vcf), collapse="\t"))
  }

  fileConn <- file(out.file)
  writeLines(output.header, fileConn)
  close(fileConn)

  print("vcf file initialized")

  fwrite(x=vcf, file=out.file, sep="\t", append = TRUE, col.names = FALSE)

  print("done")
}

