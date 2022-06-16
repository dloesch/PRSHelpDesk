#' Scale/transform PRS
#'
#' Adjusts PRS by PCA or ancestry proportions (adjusts for N-1; i.e., if 4 PCs are provided, adjusts for 3). In the case of ancestry
#' proportions, this prevents variables from summing to 1.
#' Also scales PRS within populations. If populations are not definied, populations will be determined via agglomerative clustering.
#' Note that it is assumed that the IDs for the PRS and anc.mat are matched
#' @importFrom stats lm cutree as.formula
#' @importFrom cluster agnes
#' @param IDs vector of IDS. Note: PRS and anc.mat IDs should match
#' @param PRS vector of PRS to be adjusted
#' @param anc.mat Matrix of PCAs or ancestry proportions for IDs
#' @param rescale Rescale within population?Defaults to FALSE
#' @param cluster Determine population clusters from anc.mat? Defaults to FALSE
#' @param pops Population labels for Ids. Defaults to NA. If rescaling, need to provide one of pops or cluster
#' @param k Number of clusters. Defaults to 3.
#' @export
adjustPRS <- function(IDs, PRS, anc.mat, rescale=FALSE, cluster=FALSE, pops=NULL, k=3){

  dat <- data.frame(ID=IDs, PRS=PRS, POP=pops)
  if(is.null(colnames(anc.mat))) colnames(anc.mat) <- paste0("X", 1:ncol(anc.mat))
  dat <- cbind(dat, anc.mat)

  f=as.formula(paste("PRS~", paste(colnames(dat)[4:(ncol(dat)-1)], collapse = "+")))
  fit <- lm(f, data=dat)
  PRS_adj <- fit$residuals

  #now for rescaling
  PRS_rescaled <- rep(NA, times=nrow(dat))
  if(rescale == TRUE){
    if(is.null(pops) & cluster == FALSE) stop("Either provide population labels or set cluster == TRUE")

    if(length(unique(pops)) > 1){
      for(p in unique(pops)){
        PRS_rescaled <- ifelse(dat$POP == p, scale(dat$PRS[dat$POP == p]), PRS_rescaled)
      }
    }else{
      anc.mat <- scale(anc.mat)

      print("clustering using ward method")
      hc <- cluster::agnes(anc.mat, method = "ward")

      print(paste("agglomerative coefficient:", hc$ac))
      dat$POP <- as.factor(cutree(hc, k = k))
      for(p in unique(dat$POP)){
        PRS_rescaled <- ifelse(dat$POP == p, scale(dat$PRS[dat$POP == p]), PRS_rescaled)
      }
    }
  }

  #write out
  final <- data.frame(ID=IDs, PRS=PRS, PRS_ADJ=PRS_adj, PRS_RESCALED=PRS_rescaled, stringsAsFactors = FALSE)
  return(final)
}
