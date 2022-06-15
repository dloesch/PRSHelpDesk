#' check for allele frequencies
#'
#' check if allele frequencies are already calculated
#'
#' @param ref.file file path of genotype data
#' @param bcftools path to bcftools
#' @export
checkAF <- function(ref.file, bcftools = "bcftools"){
  if(file.type(ref.file) == "bed"){
    test <- FALSE
  }else{
    p <- pipe(paste(bcftools,  "view -h", ref.file, "| grep AF"), "rt")
    test <- readLines(p)
    close(p)
    if(length(test) >0){
      test <- TRUE
    }else{
      test <- FALSE
    }
  }
  return(test)
}
