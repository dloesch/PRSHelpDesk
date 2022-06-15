#' Predict with PRS
#'
#' Generates predictions using full (covariates+PRS) and base (covariates only) models across k folds where (N-N/k) is used to predict N/k subjects
#' Predictions can then be used for generating figures, if desired (i.e. ROC)
#' @importFrom stats as.formula lm glm predict
#' @param IDs vector of IDS for PRS, trait, folds, cov matrix, etc.
#' @param PRS vector of PRS estimates
#' @param trait vector of trait. If binary, cases=1, controls=0
#' @param folds vector of assigned folds for each subject.
#' @param cov covariate matrix
#' @param binary is the trait binary? defaults to FALSE
#' @importFrom stats cor lm glm sd as.formula
#' @export
#'
predictPRS <- function(IDs, PRS, trait, folds, cov, binary=FALSE){
  dat <- data.frame(ID=IDs, PRS=PRS, trait=trait, folds=folds, stringsAsFactors = FALSE)

  dat <- cbind(dat, cov)

  f=as.formula(paste0("trait~PRS", "+", paste(colnames(dat)[5:ncol(dat)], collapse = "+")))
  b=as.formula(paste0("trait~", "+", paste(colnames(dat)[5:ncol(dat)], collapse = "+")))

  k <- max(folds)

  for(i in 1:k){
    train <- dat[dat$folds != i,]
    test <- dat[dat$folds == i,]

    if(binary == FALSE){
      new_fit <- lm(f,data=train)
      new_base <- lm(b, data=train)
    }else{
      new_fit <- glm(f,data=train, family="binomial")
      new_base <- glm(b, data=train, family="binomial")
    }

    pred <- predict(new_fit,newdata=test, type="response")
    pred2 <- predict(new_base,newdata=test, type="response")

    foo <- data.frame(ID=test$ID, TRAIT=test$trait, BASE_PRED=pred2, FULL_PRED=pred, stringsAsFactors = FALSE)
    if(i == 1){
      p <- foo
    }else{
      p <- rbind(p, foo)
    }

  }

  return(p)
}
