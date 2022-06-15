#' calculate allele counts
#'
#' Calculates allele counts of genotypes in a data.frame, matrix, or vector
#' @param geno data.frame, matrix, or vector of genotypes
#' @return vector of allele counts
#' @export
calcAC <- function(geno){
  if(is.matrix(geno) | is.data.frame(geno)){
    geno <- to_dosage(geno)
    a_counts <- rowSums(geno,na.rm=TRUE)

  }else{
    geno<- ifelse(geno== "0/0" | geno== "0|0", 0,
                  ifelse(geno== "0/1" | geno== "0|1", 1,
                         ifelse(geno== "1/0" | geno== "1|0", 1,
                                ifelse(geno== "1/1" | geno== "1|1", 2, NA))))
    a_counts <- sum(geno, na.rm = TRUE)
  }

  return(a_counts)
}
