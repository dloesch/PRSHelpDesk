#' empirical p-value for PRS
#'
#' Generates empirical p-value via N permuations. Adds 1 to both numerator and denominator to preven emp.p = 0.
#' @param PRS vector of PRS estimates
#' @param trait vector of trait
#' @param cov covariate matrix
#' @param perms Number of permuations. Defaults to 1000
#' @importFrom stats lm as.formula
#' @export
#'
empPVAL <- function(PRS, trait, cov, perms=1000){
  dat <- data.frame(PRS=PRS, trait=trait, stringsAsFactors = FALSE)

  dat <- cbind(dat, cov)

  f=as.formula(paste0("trait~PRS", "+", paste(colnames(dat)[4:ncol(dat)], collapse = "+")))

  fit <- lm(f, data = dat)

  pval <- summary(fit)$coefficients[2,4]

  ind.var <- c()
  for(i in 1:perms){
    dat$trait <- sample(dat$trait, replace = FALSE, size = nrow(dat))
    fit <- lm(f, data = dat)
    new.p <- summary(fit)$coefficients[2,4]
    ind.var <- c(ind.var, ifelse(new.p < pval, 1, 0))

  }

  emp.p <- (sum(ind.var)+1)/(perms+1)
  return(emp.p)
}
