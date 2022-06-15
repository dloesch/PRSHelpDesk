#' Match betas with target
#'
#' Match direction of effect with betas. Set betas of failed snps to 0 (i.e. not in reference, duplicates, mismatch).
#' @param mgData mgData class object
#' @param CGAT Keep CG/AT sites? Defaults to FALSE
#' @param AF Run allele frequency check? Defaults to FALSE
#' @param af.cutoff allele frequency difference cutoff. Defaults to NA (which uses 2*SD)
#' @export

matchBETA <- function(mgData, CGAT=FALSE, AF=FALSE, af.cutoff=NA){
  x <- as.data.frame(mgData)

  x$BETA <- ifelse(x$FILTER == "FAIL", 0, x$BETA)

  if(CGAT == FALSE){
    x$BETA <- ifelse(x$CGAT == TRUE, 0, x$BETA)
  }
  if(AF == TRUE){
    af.filter <- testAF(x$AF1, x$AF2, cutoff = af.cutoff, align = x$info)
    x$BETA[af.filter == "PASS"] <- 0
  }

  #adjust for allele switches
  betas <- ifelse(x$INFO == "SWITCH" | x$INFO == "FLIP_SWITCH", x$BETA*-1, x$BETA)

  return(betas)

}

