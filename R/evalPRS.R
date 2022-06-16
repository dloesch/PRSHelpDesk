#' evaluate PRS
#'
#' Performs more detailed evaluation PRS than the testPRS() function. Will report mean pvals, R2/pseudoR2, MSE/RMSE OR
#' AUC/accuracy/sensitivty/specificity (binary traits). R2/pseudoR2 is partial R2 obtained by Full.R2 - Base.R2. Pval is p-value of PRS in full model
#' @param PRS vector of PRS estimates
#' @param trait vector of trait. If binary, set cases=1, controls=0
#' @param folds vector of assigned folds for each subject.
#' @param cov covariate matrix
#' @param binary false
#' @importFrom stats lm glm predict sd formula
#' @importFrom DescTools PseudoR2
#' @importFrom pROC roc
#' @export
#'
evalPRS <- function(PRS, trait, folds, cov, binary=FALSE){
  dat <- data.frame(PRS=PRS, trait=trait, folds=folds, stringsAsFactors = FALSE)
  if(is.null(colnames(cov))) colnames(cov) <- paste0("X", 1:ncol(cov))
  dat <- cbind(dat, cov)

  f=as.formula(paste0("trait~PRS", "+", paste(colnames(dat)[4:ncol(dat)], collapse = "+")))
  b=as.formula(paste0("trait~", "+", paste(colnames(dat)[4:ncol(dat)], collapse = "+")))
  k <- max(folds)

  b.mse <- c()
  f.mse <- c()
  r2 <- c()
  pval <- c()
  f.auc <- c()
  f.sens <- c()
  f.spec <- c()
  f.bal <- c()
  f.acc <- c()
  b.auc <- c()
  b.sens <- c()
  b.spec <- c()
  b.bal <- c()
  b.acc <- c()

  for(i in 1:k){
    train <- dat[dat$folds != i,]
    test <- dat[dat$folds == i,]

    if(binary == FALSE){
      new_fit <- lm(f,data=train)
      new_base <- lm(b, data=train)

      pred <- predict(new_fit,newdata=test, type="response")
      pred2 <- predict(new_base,newdata=test, type="response")

      #MSE
      f.mse <- c(f.mse, sum((test$trait-pred))^2/nrow(test))
      b.mse <- c(b.mse, sum((test$trait-pred2)^2)/nrow(test))

      #PVAL
      pval <- c(pval, summary(new_fit)$coefficients[2,4])

      #R2
      r2 <- c(r2, summary(new_fit)$r.squared - summary(new_base)$r.squared)

    }else{
      new_fit <- glm(f,data=train, family="binomial")
      new_base <- glm(b, data=train, family="binomial")

      pred <- predict(new_fit,newdata=test, type="response")
      pred2 <- predict(new_base,newdata=test, type="response")

      #AUC
      f.auc <- c(f.auc, as.numeric(roc(test$trait, pred)$auc))
      b.auc <- c(b.auc, as.numeric(roc(test$trait, pred2)$auc))

      #Sensitivity/Specificity/Accuracy/Balanced Accuracy
      #full model
      pred <- ifelse(pred > 0.5, 1, 0)
      tp <- length(test$trait[test$trait == 1 & pred == 1])
      tn <- length(test$trait[test$trait == 0 & pred == 0])
      fp <- length(test$trait[test$trait == 0 & pred == 1])
      fn <- length(test$trait[test$trait == 1 & pred == 0])
      f.sens <- c(f.sens, tp/(tp+fn))
      f.spec <- c(f.spec, tn/(tn+fp))
      f.acc <- c(f.acc, (tp+tn)/nrow(test))
      f.bal <- c(f.bal, ((tp/(tp+fn)) + (tn/(tn+fp)))/2 )
      #base
      pred2 <- ifelse(pred2 > 0.5, 1, 0)
      tp <- length(test$trait[test$trait == 1 & pred2 == 1])
      tn <- length(test$trait[test$trait == 0 & pred2 == 0])
      fp <- length(test$trait[test$trait == 0 & pred2 == 1])
      fn <- length(test$trait[test$trait == 1 & pred2 == 0])
      b.sens <- c(b.sens, tp/(tp+fn))
      b.spec <- c(b.spec, tn/(tn+fp))
      b.acc <- c(b.acc, (tp+tn)/nrow(test))
      b.bal <- c(b.bal, ((tp/(tp+fn)) + (tn/(tn+fp)))/2 )

      #PVAL
      pval <- c(pval, summary(new_fit)$coefficients[2,4])

      #PSEUD0R2
      r2 <- c(r2, PseudoR2(new_fit, which="Nagelkerke") - PseudoR2(new_base, which="Nagelkerke"))

    }

  }
  pval.sd <- sd(pval)
  pval <- mean(pval)

  r2.sd <- sd(r2)
  r2 <- mean(r2)

  if(binary == TRUE){
    #get mean and sd of each metric
    f.auc.sd <- sd(f.auc)
    f.auc <- mean(f.auc)
    f.sens.sd <- sd(f.sens)
    f.sens <- mean(f.sens)
    f.spec.sd <- sd(f.spec)
    f.spec <- mean(f.spec)
    f.acc.sd <- sd(f.acc)
    f.acc <- mean(f.acc)
    f.bal.sd <- sd(f.bal)
    f.bal <- mean(f.bal)
    b.auc.sd <- sd(b.auc)
    b.auc <- mean(b.auc)
    b.sens.sd <- sd(b.sens)
    b.sens <- mean(b.sens)
    b.spec.sd <- sd(b.spec)
    b.spec <- mean(b.spec)
    b.acc.sd <- sd(b.acc)
    b.acc <- mean(b.acc)
    b.bal.sd <- sd(b.bal)
    b.bal <- mean(b.bal)

    #final results
    results <- data.frame(MODEL=c("BASE", "FULL"),
                          R2=c(NA, r2), R2.SD=c(NA, r2.sd), PVAL=c(NA, pval), PVAL.SD=c(NA, pval.sd),
                          ACC=c(b.acc, f.acc), ACC.SD=c(b.acc.sd, f.acc.sd),
                          BAL=c(b.bal, f.bal), BAL.SD=c(b.bal.sd, f.bal.sd),
                          SENS=c(b.sens, f.sens), SENS.SD=c(b.sens.sd, f.sens.sd),
                          SPEC=c(b.spec, f.spec), SPEC.SD=c(b.spec.sd, f.spec.sd),
                          stringsAsFactors = FALSE)
  }else{
    #rmse
    b.rmse <- sqrt(b.mse)
    f.rmse <- sqrt(f.mse)

    #get mean and sd of each metric
    b.mse.sd <- sd(b.mse)
    b.mse <- mean(b.mse)
    f.mse.sd <- sd(f.mse)
    f.mse <- mean(f.mse)
    b.rmse.sd <- sd(b.rmse)
    b.rmse <- mean(b.rmse)
    f.rmse.sd <- sd(f.rmse)
    f.rmse <- mean(f.rmse)

    results <- data.frame(MODEL=c("BASE", "FULL"),
                          R2=c(NA, r2), R2.SD=c(NA, r2.sd), PVAL=c(NA, pval), PVAL.SD=c(NA, pval.sd),
                          MSE=c(b.mse, f.mse), MSE.SD=c(b.mse.sd, f.mse.sd),
                          RMSE=c(b.rmse, f.rmse), RMSE.SD=c(b.rmse.sd, f.rmse.sd),
                          stringsAsFactors = FALSE)
  }

  return(results)
}
