#' get dosages from imputed data
#'
#' If data is imputed, dosages will be included
#' @param geno matrix or data frame of imputed genotypes
#' @return matrix of imputed genotypes
#' @export
imputed_dosage <- function(geno){
  if(!is.matrix(geno)){
    geno <- as.matrix(geno)
  }
  geno <- apply(X = geno, c(1,2), FUN = function(x) as.numeric(unlist(unname(strsplit(x, split = ":")))[2]))

  return(geno)
}


