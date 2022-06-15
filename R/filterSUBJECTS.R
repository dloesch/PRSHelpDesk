#' Filter subjects
#'
#' Filter subjects by providing input data for each subject (i.e. heterozygosity, PC, ancestry proportions, etc.). Subjects will be idenfited that
#' exceed 3*SD from the mean of the provided data. Note that filterRELS() identifies relatives
#' @importFrom stats sd
#' @param IDs vector of ids
#' @param INPUT vector of values for filtering
#' @export

filterSUBJECTS <- function(IDs, INPUT){
  dat <- data.frame(ID=IDs, y=INPUT)

  upper <- mean(dat$y, na.rm = TRUE) + 3*sd(dat$y, na.rm = TRUE)
  lower <- mean(dat$y, na.rm = TRUE) - 3*sd(dat$y, na.rm = TRUE)

  dat <- dat[dat$y < lower | dat$y > upper,]

  return(dat$ID)
}

