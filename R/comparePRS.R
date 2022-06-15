#' Compare PRS models
#'
#' Compares PRS models using the Williams test from the psych package as described by Pain et al. 2021. Will report Base vs PRS1,
#' Base vs. PRS2, PRS1 vs. PRS2
#' @importFrom stats as.formula
#' @param PRS1 first PRS to be tested
#' @param PRS2 second PRS to be tested
#' @param trait vector of trait. If binary, cases=1, controls=0
#' @param folds vector of assigned folds for each subject.
#' @param cov covariate matrix
#' @param binary is the trait binary? defaults to FALSE
#' @importFrom stats cor lm glm sd
#' @importFrom psych paired.r
#' @export
#'
comparePRS <- function(PRS1, PRS2, trait, folds, cov, binary=FALSE){
  dat <- data.frame(PRS1=PRS1, PRS2=PRS2, trait=trait, folds=folds, stringsAsFactors = FALSE)

  dat <- cbind(dat, cov)

  f1=as.formula(paste0("trait~PRS1", "+", paste(colnames(dat)[5:ncol(dat)], collapse = "+")))
  f2=as.formula(paste0("trait~PRS2", "+", paste(colnames(dat)[5:ncol(dat)], collapse = "+")))
  b=as.formula(paste0("trait~", "+", paste(colnames(dat)[5:ncol(dat)], collapse = "+")))

  k <- max(folds)
  f1.cor <- c()
  f2.cor <- c()
  b.cor <- c()
  prs1_prs2 <- c()
  prs1_base <- c()
  prs2_base <- c()
  for(i in 1:k){
    train <- dat[dat$folds != i,]
    test <- dat[dat$folds == i,]

    if(binary == FALSE){
      prs1_fit <- lm(f1,data=train)
      prs2_fit <- lm(f2,data=train)
      base <- lm(b, data=train)
    }else{
      prs1_fit <- glm(f1,data=train, family="binomial")
      prs2_fit <- glm(f2,data=train, family="binomial")
      base <- glm(b, data=train, family="binomial")
    }

    pred1 <- predict(prs1_fit,newdata=test, type="response")
    pred2 <- predict(prs2_fit,newdata=test, type="response")
    pred3 <- predict(base,newdata=test, type="response")

    prs1_base <- c(prs1_base, cor(pred1, pred3))
    prs2_base <- c(prs2_base, cor(pred2, pred3))
    prs1_prs2 <- c(prs1_prs2, cor(pred1, pred2))

    f1.cor <- c(f1.cor, cor(pred1, test$trait))
    f2.cor <- c(f2.cor, cor(pred2, test$trait))
    b.cor <- c(b.cor, cor(pred3, test$trait))

  }
  f1.sd <- sd(f1.cor)
  f1.cor <- mean(f1.cor)
  f2.sd <- sd(f2.cor)
  f2.cor <- mean(f2.cor)
  base.sd <- sd(b.cor)
  b.cor <- mean(b.cor)

  #standard error
  N <- length(trait)
  f1.se <- (1-f1.cor^2)/sqrt(N-2)
  f2.se <- (1-f1.cor^2)/sqrt(N-2)
  base.se <- (1-b.cor^2)/sqrt(N-2)

  #get pvals
  prs1_base <- mean(prs1_base)
  prs2_base <- mean(prs2_base)
  prs1_prs2 <- mean(prs1_prs2)

  prs1_base.pval <- paired.r(f1.cor, b.cor, prs1_base, n=N, twotailed = TRUE)$p
  prs2_base.pval <- paired.r(f2.cor, b.cor, prs2_base, n=N, twotailed = TRUE)$p
  prs1_prs2.pval <- paired.r(f1.cor, f2.cor, prs1_prs2, n=N, twotailed = TRUE)$p

  #correlation between PRS
  PRS1_PRS2.cor <- cor(PRS1,PRS2)

  final <- data.frame(N=N,
                      R.PRS1=f1.cor, SD.PRS1= f1.sd, SE.PRS1=f1.se,
                      R.PRS2=f2.cor, SD.PRS2= f2.sd, SE.PRS2=f2.se,
                      R.BASE=b.cor, SD.BASE= base.sd, SE.BASE=base.se,
                      P.BASE_PRS1= prs1_base.pval,
                      P.BASE_PRS2= prs2_base.pval,
                      P.PRS1_PRS2= prs1_prs2.pval,
                      R.PRS1_PRS2= PRS1_PRS2.cor,
                      stringsAsFactors = FALSE)

  return(final)

}

