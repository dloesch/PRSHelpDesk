#' Flip dosages
#'
#' "Flips" dosages when allele switch occurs
#' @param dosage vector or matrix of genotype dosages
#' @export
flip_dosage <- function(dosage){
  dosage <- ifelse(dosage ==2, 0,
                   ifelse(dosage ==1, 1,
                   ifelse(dosage == 0, 2, NA)))
  return(dosage)
}
