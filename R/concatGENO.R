#' Concatenate genotype files
#'
#' Wrapper function for concatenating genotype files
#' @param geno.list Text file with genotype filepaths or vector with file names
#' @param prefix Prefix, defaults to output
#' @param file.type Select one of vcf, bed, bcf
#' @param plink Path to plink
#' @param bcftools Path to bcftools
#' @importFrom utils write.table
#' @export

concatGeno <- function(geno.list, prefix="output", file.type=c("vcf", "bed", "bcf"),
                       plink="plink", bcftools="bcftools"){
  #run some checks
  if(length(geno.list) ==1){
    if(!file.exists(geno.list)) stop("Unable to read file")
  }else{
    if(!is.vector(geno.list)) stop("Provide character vector of filepaths")
    for(f in geno.list){
      if(!file.exists(f)) stop("Unable to read genotype files")
    }
    temp <- tempfile()
    write.table(geno.list, temp, col.names = FALSE, row.names = FALSE, sep = '\t', quote=FALSE)
    geno.list <- temp
  }

  if(length(file.type) > 1){
    stop("Please select one file type, vcf, bed, or bcf")
  }else{
    file.type <- tolower(file.type)
  }

  if(file.type == "bed"){
    system2(plink, paste("--merge-list", geno.list, "--make-bed --out", prefix))
  }else if(file.type == "vcf"){
    out.file <- paste0(prefix, ".vcf.gz")
    system2(bcftools, paste("concat -f", geno.list, "-Oz -o", out.file))
  }else if(file.type == "bcf"){
    out.file <- paste0(prefix, ".bcf")
    system2(bcftools, paste("concat -f", geno.list, "-Ob -o", out.file))
  }else{
    stop("Genotype files need to be bed, bcf, or vcf")
  }
}
