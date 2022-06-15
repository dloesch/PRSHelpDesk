#' calculate number of chromosomal observations
#'
#' calculate number of chromosomal observations in a data.frame, matrix, or vector
#' @param geno data.frame, matrix, or vector of genotypes
#' @return vector of chromosomal observations
#' @export
calcNobs <- function(geno){
  geno <- to_dosage(geno)
  geno[!is.na(geno)] <- 1
  a_counts <- rowSums(geno,na.rm=TRUE)
  Nobs <- rowSums(geno, na.rm=TRUE)*2
  return(Nobs)
}
