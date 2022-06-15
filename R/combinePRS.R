#' Combine PRS
#'
#' Combines PRS using provided weights/coefficients and intercept
#'
#' @param IDs vector of IDS matching PRS1,PRS2, etc.
#' @param PRS1 PRS1 vector
#' @param w1 weight/coefficient for PRS1
#' @param PRS2 PRS2 vector
#' @param w2 weight/coefficient for PRS2
#' @param PRS3 PRS3 vector (optional)
#' @param w3 weight/coefficient for PRS3 (optional)
#' @param PRS4 PRS4 vector (optional)
#' @param w4 weight/coefficient for PRS4 (optional)
#' @param intercept intercept for linear combination. Defaults to 0
#' @export

combinePRS <- function(IDs, PRS1, w1, PRS2, w2, PRS3=0, w3=0,PRS4=0, w4=0, intercept=0){
  PRS_COMB <- PRS1*w1 + PRS2*w2 + PRS3*w3 + PRS4*w4 + intercept
  dat <- data.frame(ID=IDs, PRS_COMB=PRS_COMB, stringsAsFactors = FALSE)

  return(dat)
}
