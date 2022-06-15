#' calculate allele frequencies
#'
#' Calculates allele frequencies of genotypes in a data.frame, matrix, or vector
#' @param geno data.frame, matrix, or vector of genotypes
#' @return vector of allele frequencies
#' @export
calcAF <- function(geno){
  if(is.matrix(geno) | is.data.frame(geno)){
    geno <- to_dosage(geno)
    #second matrix to just count number of observations in each row
    geno2 <- geno
    geno2[!is.na(geno2)] <- 1
    af_counts <- rowSums(geno,na.rm=TRUE)
    Nobs <- rowSums(geno2, na.rm=TRUE)*2
    af <- af_counts/Nobs
  }else{
    geno<- ifelse(geno== "0/0" | geno== "0|0", 0,
                  ifelse(geno== "0/1" | geno== "0|1", 1,
                         ifelse(geno== "1/0" | geno== "1|0", 1,
                                ifelse(geno== "1/1" | geno== "1|1", 2, NA))))
    af <- sum(geno, na.rm = TRUE)/(length(geno[!is.na(geno)])*2)
  }

  return(unlist(unname(af)))
}
