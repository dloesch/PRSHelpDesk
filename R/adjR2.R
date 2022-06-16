#' adjusts R2
#'
#' For binary traits, estimates R2 on the liability scale, pseudoR2 on the liability scale (and probit/logit) from Lee et al. 2012 "Better Coefficient".
#' For both quantitative and binary traits, if PRS is the result of a linear combination of multiple PRS, computes adjusted R2 accounting for greater model complexity
#' (Marquez-Luna et al. 2017)
#' @param PRS vector of PRS estimates
#' @param trait vector of trait. If binary, Cases=1, controls=0
#' @param cov covariate matrix
#' @param nPRS number of PRS in combination
#' @param prev trait prevalence
#' @param binary Binary or quantitative trait. Defaults to TRUE
#' @importFrom stats glm lm as.formula qnorm dnorm var binomial
#' @export
#'
adjR2 <- function(PRS, trait, cov, nPRS=1, prev=NA, binary=TRUE){
  dat <- data.frame(PRS=PRS, trait=trait, stringsAsFactors = FALSE)
  if(is.null(colnames(cov))) colnames(cov) <- paste0("X", 1:ncol(cov))
  dat <- cbind(dat, cov)

  if(binary == TRUE){
    f=as.formula(paste0("trait~PRS", "+", paste(colnames(dat)[4:ncol(dat)], collapse = "+")))
    b=as.formula(paste0("trait~", "+", paste(colnames(dat)[4:ncol(dat)], collapse = "+")))

    full <- glm(f, data = dat, family = "binomial")
    base <- glm(b, data = dat, family = "binomial")

    nt <- nrow(dat)
    ncase <- nrow(dat[dat$trait == 1,])
    ncont <- nrow(dat[dat$trait == 0,])
    if(is.na(prev)) stop("specify prevalence rate")
    K <- prev
    P <- ncase/nt
    thd <- qnorm(1-K)
    zv <- dnorm(thd)
    mv <- zv/K
    mv2 <- -mv*K/(1-K)

    R2O <- var(full$fitted.values)/(ncase/nt*ncont/nt) - var(base$fitted.values)/(ncase/nt*ncont/nt)
    theta=mv*(P-K)/(1-K)*(mv*(P-K)/(1-K) -thd)
    cv=K*(1-K)/zv^2*K*(1-K)/(P*(1-P))
    R2=R2O*cv/(1+R2O*theta*cv)

    #logistic liability scale
    full <- glm(f, family=binomial(logit), data = dat)
    base <- glm(b, family=binomial(logit), data = dat)
    r2.logit <- var(full$linear.predictors)/(var(full$linear.predictors)+pi^2/3) - var(base$linear.predictors)/(var(base$linear.predictors)+pi^2/3)

    #probit liability scale
    full <- glm(f, family=binomial(probit), data = dat)
    base <- glm(b, family=binomial(probit), data = dat)
    r2.probit <- var(full$linear.predictors)/(var(full$linear.predictors)+1) - var(base$linear.predictors)/(var(base$linear.predictors)+1)

    results <- data.frame(R2.liab= R2, pR2.logit=r2.logit, pR2.probit=r2.probit, stringsAsFactors = FALSE)

  }else{
    full <- lm(f, data = dat)
    base <- lm(b, data = dat)
    R2 <- summary(full)$r.squared - summary(base)$r.squared
    results <- data.frame(R2=R2, stringsAsFactors = FALSE)
  }

  if(nPRS > 1){
    p <- nPRS
    N <- nrow(dat)
    R_adj <- R2 - (1-R2)*(p/(N-p-1))
    results$R_ADJ <- R_adj
  }

    return(results)
}
