#' flip allele to opposite strand
#'
#'@param x A,T,C, or G
#'@return flipped allele
#'@examples
#' flip_allele("A")
#'@export
flip_allele <- function(x){
  x <- toupper(x)
  y <- ifelse(x == "A", "T", ifelse(x == "T", "A", ifelse(x == "C", "G", ifelse(x == "G", "C", NA))))
  return(y)
}
