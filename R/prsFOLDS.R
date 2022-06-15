#' PRS folds
#'
#' Generates folds for PRS folds
#'
#' For quantitative traits, assigns subjects to k folds randomly. For binary traits, assigns equal number of case/controls to k folds.
#' @param subjects vector of subject IDS
#' @param binary Binary trait for stratified sampling. Defaults to FALSE
#' @param case_control Vector of case/control status. Defaults to NA
#' @param folds Number of folds. Defaults to 10
#' @export
prsFOLDS <- function(subjects,binary=FALSE, case_control=NA, folds=10){
  data <- data.frame(ID=sample(subjects, length(subjects)), stringsAsFactors = FALSE)
  k <- folds
  if(binary == FALSE){
    data$FOLD <- rep(1:k)
  }else{
    if(is.na(case_control)) stop("provide vector of case/control membership")
    data$STATUS <- case_control
    groups <- unique(case_control)
    #split
    group_1 <- data$ID[data$STATUS == groups[1]]
    group_2 <- data$ID[data$STATUS == groups[2]]
    N1 <- round(length(group_1)/k)
    N2 <- round(length(group_2)/k)

    #assign
    data$FOLD <- k
    for(i in 1:(k-1)){
      fold <- c(sample(x=group_1, size = N1, replace = FALSE),
                sample(x=group_2, size = N2, replace = FALSE))
      data$FOLD <- ifelse(data$ID %in% fold, i, data$FOLD)
      group_1 <- group_1[!group_1 %in% fold]
      group_2 <- group_2[!group_2 %in% fold]
    }
  }

}
