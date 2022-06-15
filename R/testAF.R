#' Test allele frequencies
#'
#' Simple test of allele frequencies. Returns a vector indicating PASS/FAIL.
#' SNP is flagged as FAIL if abs(AF DIFF) > mean+2*SD. Only variants with mean difference < 0.5 are included in estimation of mean
#' If desired, user can provide a hard cutoff, i.e. an allele frequency difference of 0.2 (absolute value)
#' @importFrom ggplot2 ggplot aes geom_point theme_bw ylab xlab theme
#' @importFrom MASS kde2d
#' @importFrom stats density
#' @importFrom rlang .data
#' @param af1 alelle frequencies from target cohort
#' @param af2 allele frequencies from GWAS cohort
#' @param cutoff Desired abs(allele frequency difference) cutoff
#' @param align vector indicating SAME/SWITCH/FLIP/FLIP_SWITCH. Defaults to NA
#' @param plot Generates plot of allele frequencies. Defaults to FALSE
#' @export

testAF <- function(af1, af2, cutoff=NA, align=NA, plot=FALSE){
  if(length(align) == length(af1)){
    af1 <- ifelse(align == "SWITCH" | align == "FLIP_SWITCH", 1 -af1, af1)
  }

  diff <- abs(af1 - af2)
  if(is.na(cutoff)) cutoff <- mean(diff[diff < 0.5]) + 2*sd(diff[diff < 0.5])
  cutoff <- as.numeric(cutoff)
  if(cutoff <=0 | cutoff >=1) stop("Cutoff needs to be >0 & <1")
  af_filter <- ifelse(diff > cutoff, "FAIL", "PASS")

  if(plot == TRUE){
    #code from https://slowkow.com/notes/ggplot2-color-by-density/
    get_density <- function(x, y, ...) {
      dens <- kde2d(x, y, ...)
      ix <- findInterval(x, dens$x)
      iy <- findInterval(y, dens$y)
      ii <- cbind(ix, iy)
      return(dens$z[ii])
    }
    temp$density <- get_density(temp$AF1, temp$AF2, n=100)

    g <- ggplot(data=temp, aes(.data$AF1, .data$AF2, color=density)) +
      geom_point(size=1)+
      theme_bw()
    g <-  g + theme(legend.position='none')
    g <- g+ xlab("TARGET COHORT AF") + ylab("GWAS COHORT AF")
    print(g)
  }

  return(af_filter)

}
