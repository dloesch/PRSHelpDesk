#' Convert genotypes to dosages
#'
#' Converts genotypes, i.e. 1/0, to genotype dosages, i.e 1.
#' @param geno Genotype data in a matrix, data.frame, or vector
#' @return vector or matrix with genotype dosages
#' @export
to_dosage <- function(geno){
  #convert genotypes to genotype dosages
  if(is.data.frame(geno)){
    geno <- as.matrix(geno)
  }

  #if imputed, get hard-coded genotype
  if(length(geno[grep(":",geno)]) > 0){geno <- substr(geno, 1, regexpr("\\:", geno)-1)}
  geno <- ifelse(geno == "1|1" | geno == "1/1", 2,
                 ifelse(geno == "1|0" | geno == "1/0" | geno == "0/1" | geno == "0/1", 1,
                        ifelse(geno == "0|0" | geno == "0/0", 0, NA)))


  return(geno)
}
