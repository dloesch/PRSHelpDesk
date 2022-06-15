#' Learn weights for linear combination of PRS
#'
#' Learn weights for linear combination of PRS. Two methods will be implemented: 1) fitting a regression model and 2)
#' learning the weights empirically where the weights sum to 1. The latter method appears to work better for admixed cohorts.
#'
#' @importFrom data.table setDT
#' @importFrom stats lm glm cor as.formula
#' @param PRS1 vector for PRS1
#' @param PRS2 vector for PRS2
#' @param PRS3 vector for PRS3 (optional)
#' @param PRS4 vector for PRS4 (optional)
#' @param trait vector of trait. If binary, cases=1, controls=0
#' @param folds vector of assigned folds for each subject.
#' @param cov covariate matrix
#' @param binary is the trait binary? defaults to FALSE
#' @export
#'
weightPRS <- function(PRS1, PRS2, PRS3=NA, PRS4=NA, trait, folds, cov, binary=FALSE){
  W <- X <- Y <- Z <- NULL
  dat <- data.frame(PRS1=PRS1, PRS2=PRS2, PRS3=PRS3, PRS4=PRS4, trait=trait, folds=folds, stringsAsFactors = FALSE)

  dat <- cbind(dat, cov)

  if(is.na(PRS3) & is.na(PRS4)){
    #obtain coefficients from Y~PRS1+PRS2
    fit <- lm(trait~PRS1+PRS2)
    intercept <- fit$coefficients[1]
    w1 <- fit$coefficients[2]
    w2 <- fit$coefficients[3]
    w3 <- NA
    w4 <- NA
    #obtain weights empirically
    f=as.formula(paste0("trait~PRS_SUM", "+", paste(colnames(dat)[7:ncol(dat)], collapse = "+")))
    dat$PRS3 <- 0
    dat$PRS4 <- 0
    #set up datatable of possible weights to iterate over
    df <- data.frame(W=(1:100)/100)
    setDT(df)
    df <- df[W >= 0.1 ]
    df$X <- rep(1, times=nrow(df)) - df$W
    df <- df[X >= 0.1]
    df$Y <- 0
    df$Z <- 0

  }else if(!is.na(PRS3) & is.na(PRS4)){
    #obtain coefficients from Y~PRS1+PRS2+PRS3
    fit <- lm(trait~PRS1+PRS2+PRS3)
    intercept <- fit$coefficients[1]
    w1 <- fit$coefficientsp[2]
    w2 <- fit$coefficients[3]
    w3 <- fit$coefficients[4]
    w4 <- NA
    #obtain weights empirically
    f=as.formula(paste0("trait~PRS_SUM", "+", paste(colnames(dat)[7:ncol(dat)], collapse = "+")))


    #set up datatable of possible weights to iterate over
    df <- expand.grid(W=(1:100)/100, X=(1:100)/100)
    setDT(df)
    df <- df[(W+X) < 1]
    df$Y <- rep(1, times=nrow(df)) - (df$W + df$X)
    df <- df[W >= 0.1 ]
    df <- df[X >= 0.1]
    df <- df[Y >= 0.1]
    df$Z <- 0
  }else{
    #obtain coefficients from Y~PRS1+PRS2+PRS3+PRS4
    fit <- lm(trait~PRS1+PRS2+PRS3+PRS4)
    intercept <- fit$coefficients[1]
    w1 <- fit$coefficientsp[2]
    w2 <- fit$coefficients[3]
    w3 <- fit$coefficients[4]
    w4 <- fit$coefficients[5]
    #obtain weights empirically
    f=as.formula(paste0("trait~PRS1+PRS2+PRS3+PRS4", "+", paste(colnames(dat)[7:ncol(dat)], collapse = "+")))

    df <- expand.grid(W=(1:100)/100, X=(1:100)/100, Y=(1:100)/100)
    setDT(df)
    df <- df[(W+X+Y) < 1]
    df <- df[W >= 0.1]
    df <- df[X >= 0.1 ]
    df <- df[Y >= 0.1]
    df$Z <- rep(1, times=nrow(df)) - (df$W + df$X + df$Y)
    df <- df[Z >= 0.1]
  }

  df <- as.data.frame(df)
  df$R_MEAN <- NA
  k <- max(folds)


  for(i in 1:nrow(df)){

    #linear combination
    dat$PRS_SUM <- dat$PRS1*df[i,1] + dat$PRS2*df[i,2] + dat$PRS3*df[i,3] + dat$PRS4*df[i,4]

    #correlations
    f.cor <- c()

    for(j in 1:k){
      train <- dat[dat$folds != j,]
      validation <- dat[dat$folds == j,]

      if(binary == TRUE){
        new_fit <- glm(f,data=train, family="binomial")
      }else{
        new_fit <- lm(f,data=train)
      }

      #predictions
      pred <- predict(new_fit,newdata=validation, type="response")


      #correlations
      f.cor <- c(f.cor, cor(validation$trait, pred))
    }

    #mean r
   df$R_MEAN[i] <- mean(f.cor)

  }

  w1.emp <- df$W[df$R_MEAN == max(df$R_MEAN)][1]
  w2.emp <- df$X[df$R_MEAN == max(df$R_MEAN)][1]
  w3.emp <- df$Y[df$R_MEAN == max(df$R_MEAN)][1]
  w4.emp <- df$Z[df$R_MEAN == max(df$R_MEAN)][1]

  weights <- data.frame(LABEL=c("INTERCEPT", "PRS1", "PRS2", "PRS3", "PRS4"),
                        WEIGHT=c(intercept, w1, w2, w3,w4),
                        EMP_WEIGHT=c(NA, w1.emp, w2.emp, w3.emp, w4.emp))
  return(weights)

}
