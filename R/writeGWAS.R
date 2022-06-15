#' export parsed GWAS sumstats
#'
#' write out parsed GWAS sumstats in formats for PRSice2, PRS-CS, and LD score regression
#' @param sumstats gData object of GWAS data
#' @param type PRS-CS, PRSice2, or LDSR. Defaults to PRS-cS
#' @param prefix Desired prefix. Defaults to "output"
#' @param by.chr Write our by chromosome? Defaults to FALSE
#' @importFrom data.table fwrite
#' @export

writeGWAS <- function(sumstats, type="PRS-CS", prefix="output", by.chr=FALSE){
  if(type == "PRS-CS"){
    data <- data.frame(SNP=sumstats$id, A1=sumstats$alt, A2=sumstats$ref,
                       BETA=sumstats$beta, P=sumstats$pval, stringsAsFactors = FALSE)
  }else if(type == "PRSice2"){
  data <- data.frame(SNP=sumstats$id, CHR=sumstats$chr, AF=sumstats$af, A1=sumstats$alt, A2=sumstats$ref,
                     BETA=sumstats$beta, SE=sumstats$se, P=sumstats$pval, stringsAsFactors = FALSE)
  data$CHR <- gsub("chr", "", data$CHR)
  }else if(type == "LDSR"){
    data <- data.frame(SNP=sumstats$id, A1=sumstats$alt, A2=sumstats$ref, Z=sumstats$beta/sumstats$se, N=sumstats$N,
                       stringsAsFactors = FALSE)
  }else{
    stop("Only supports PRS-CS, PRSice2, and LDSR formats")
  }
  if(by.chr == TRUE){
    for(chr in 1:22){
      fwrite(data[sumstats$chr == chr,], paste0(prefix, ".chr", chr, ".", type, ".txt"),
             sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
    }
  }else{
    fwrite(data, paste0(prefix, ".", type, ".txt.gz"),
           sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE, compress = "gzip")
  }

}
