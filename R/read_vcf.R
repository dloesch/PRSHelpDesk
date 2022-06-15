#' Read in VCF/BCF/BED
#'
#' Reads VCF/BCF/BED record into memory
#' @importFrom data.table fread
#' @importFrom utils read.table
#' @param input.vcf filepath to genotype data
#' @param plink path to plink. Defaults to working directory
#' @param bcftools path to bcftools Defaults to working directory
#' @return data.frame with vcf
#' @export
#' @examples
#' \dontrun{
#' read_vcf("data.vcf.gz")
#' }
read_vcf <- function(input.vcf, plink="plink", bcftools="bcftools"){
  if(!file.exists(input.vcf)) stop("error: unable to read file")
  if(file.type(input.vcf) == "vcf"){
    vcf <- fread(input.vcf, data.table = FALSE)
  }
  else if(file.type(input.vcf) == "bed"){
    #convert to VCF for reading
    input.vcf <- sub(".bed", "", input.vcf)
    tempvcf <- tempfile()
    system2(plink,
            paste("--bfile", input.vcf, "--recode vcf-iid bgz  --out", tempvcf),
            stdout = FALSE)
    vcf <- fread(paste0(tempvcf, "vcf.gz"), data.table=FALSE)
    system2("rm", paste0(tempvcf,"*"))
  }else if(file.type(input.vcf) == "bcf"){
    #pipe in bcf file
    bcf <- pipe(paste(bcftools, "view", input.vcf), "rt")
    vcf <- read.table(bcf, stringsAsFactors = FALSE)
    close(bcf)
  }else{
    stop("VCF/BCF/BED formats are supported")
  }
  return(vcf)
}
