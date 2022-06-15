#' Simulate PRS and a trait
#'
#' Simple function for generating simulated traits for testing purposes
#' @importFrom stats rnorm
#' @param N Number of subjects. Defaults to 500
#' @param binary Binary or quantitative trait. Defaults to TRUE
#' @export

simPRS <- function(N=500, binary=TRUE){

  PRS <- rnorm(N)
  SEX <- sample(1:2, replace=TRUE, size = N)
  AGE <- sample(30:80, replace = TRUE, size=N)

  #adjust to establish association
  if(binary == TRUE){
    trait <- sample(c(0,1), replace=TRUE, size = N)
    PRS <- ifelse(trait == 1, PRS+0.25, PRS-0.25)
    AGE <- ifelse(trait == 1, AGE+3, AGE-3)
  }else{
    PRS <- ifelse(trait > 1, PRS+0.5, ifelse(trait < 1, PRS-0.5, PRS))
    AGE <- ifelse(trait > 1, AGE+3, ifelse(trait < 1, AGE-3, AGE))
  }

  sim <- data.frame(ID=paste0("S", 1:N), TRAIT=trait, SEX=SEX, AGE=AGE, PRS=PRS, stringsAsFactors = FALSE)
  return(sim)
}
