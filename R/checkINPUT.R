#' Check input
#'
#' Checks phenotype data input (trait, PRS, covariate matrix, ancestry matrix). Creates and returns a dataframe.
#' @param trait vector of trait values
#' @param PRS vector of PRS values
#' @param folds vector of fold membership for cross-validation. Defaults to NULL
#' @param IDs vector of subject IDs. Defaults to NULL
#' @param cov covariate matrix. Can be matrix or data.frame. Defaults to NULL
#' @param anc.mat ancestry matrix. Can be PCA or ancestry proportions. Defaults to NULL
#' @param pops vector of population labels. Defaults to NULL
#' @export
checkINPUT <- function(trait, PRS, IDs=NULL, folds=NULL, cov=NULL, anc.mat=NULL, pops=NULL){
  #check trait
  if(!is.vector(trait)) stop("trait is not a vector")

  #check PRS
  if(!is.vector(PRS)) stop("PRS is not a vector")
  if(length(trait) != length(PRS)) stop("number of subjects do not match for trait and PRS")
  if(!is.numeric(PRS)) stop("PRS needs to be numeric")

  #check folds
  if(!is.null(folds)){
    if(!is.vector(folds)) stop("PRS is not a vector")
    if(length(trait) != length(folds)) stop("number of subjects do not match for trait and folds")
    if(!is.integer(folds)) stop("folds need to be integers between 1 and a resonable value (i.e. 10)")
    if(max(folds) > 10) warning("Using number of folds > 10")

    #make data.frame with folds
    dat <- data.frame(PRS=PRS, trait=trait, folds=folds, stringsAsFactors = FALSE)
  }else{
    #make data.frame without folds
    dat <- data.frame(PRS=PRS, trait=trait, stringsAsFactors = FALSE)
  }

  if(!is.null(IDs)){
    if(!is.vector(IDs)) stop("ID vector is not a vector")
    if(length(trait) != length(pops)) stop("number of subjects do not match for trait and ID")
    cols <- colnames(dat)
    dat$ID <- IDs
    dat <- dat[c("ID", cols)]
  }

  #check population vector
  if(!is.null(pops)){
    if(!is.vector(pops)) stop("population vector is not a vector")
    if(length(trait) != length(pops)) stop("number of subjects do not match for trait and pop labels")
    dat$POP <- pops
  }

  #check covariate matrix
  if(!is.null(cov)){
    if(!is.matrix(cov) & !is.data.frame(cov)) stop("covariates need to be provided as a matrix or data.frame")
    if(nrow(cov) != length(trait)) stop("number of subjects do not match for trait and covariates")
    if(is.null(colnames(cov))) colnames(cov) <- paste0("X", 1:ncol(cov))

    #add to data.frame
    dat <- cbind(dat, cov)
  }

  #check ancestry matrix
  if(!is.null(anc.mat)){
    if(!is.matrix(anc.mat) & !is.data.frame(anc.mat)) stop("covariates need to be provided as a matrix or data.frame")
    if(nrow(anc.mat) != length(trait)) stop("number of subjects do not match for trait and ancestry matrix")
    if(is.null(colnames(anc.mat))) colnames(anc.mat) <- paste0("X", 1:ncol(anc.mat))
    #add to data.frame
    dat <- cbind(dat, anc.mat)
  }

  return(dat)
}
