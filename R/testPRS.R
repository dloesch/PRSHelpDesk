#' test PRS
#'
#' Performs quick test of PRS by reporting mean Pearson's R between the trait and the prediction across K folds
#' @param PRS vector of PRS estimates
#' @param trait vector of trait. If binary, cases=1, controls=0
#' @param folds vector of assigned folds for each subject.
#' @param cov covariate matrix
#' @param binary is the trait binary? defaults to FALSE
#' @importFrom stats cor lm glm sd as.formula
#' @export
#'
testPRS <- function(PRS, trait, folds, cov, binary=FALSE){
  dat <- data.frame(PRS=PRS, trait=trait, folds=folds, stringsAsFactors = FALSE)

  dat <- cbind(dat, cov)

  f=as.formula(paste0("trait~PRS", "+", paste(colnames(dat)[4:ncol(dat)], collapse = "+")))
  b=as.formula(paste0("trait~", "+", paste(colnames(dat)[4:ncol(dat)], collapse = "+")))
  k <- max(folds)
  f.cor <- c()
  b.cor <- c()
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
    f.cor <- c(f.cor, cor(pred, test$trait))
    b.cor <- c(b.cor, cor(pred2, test$trait))
  }
  full.sd <- sd(f.cor)
  full <- mean(f.cor)
  base.sd <- sd(b.cor)
  base <- mean(b.cor)

  #improvement
  imp <- ((full - base)/base)*100

  #standard error
  N <- length(trait)
  full.se <- (1-full^2)/sqrt(N-2)
  base.se <- (1-base^2)/sqrt(N-2)

  final <- data.frame(N=N,R.FULL=full, SD.FULL= full.sd, SE.FULL=full.se, R.BASE=base, SD.BASE= base.sd, SE.BASE=base.se, IMP=imp,
                      stringsAsFactors = FALSE)

  return(final)
}
