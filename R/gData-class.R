#' Class gData
#'
#' Constructor function and related methods for S3 class gData
#'
#' @param chr chromosome vector
#' @param pos position vector in bp
#' @param id snp ids vector
#' @param ref ref/a2 allele vector
#' @param alt alt/a1 allele vector
#' @param N sample size vector
#' @param beta beta effect size, if data is from summary statistics. Defaults to NA
#' @param se standard error of the beta. Defaults to NA.
#' @param pval pvalue, if data is from summary statistics. Defaults to NA
#' @param af allele frequencies. Defaults to NA
#' @param build Genome build. Defaults to hg19
#' @export
gData <- function(chr, pos, id, ref, alt, N, beta=NA, se=NA, pval=NA, af=NA, build="hg19"){
  #check build
  build <- tolower(build)
  if(build != "hg19" & build != "hg38")  stop("build must be hg19 or hg38")
  #remove chr label, if needed
  chr<- gsub("chr", "", chr)
  chr <- as.numeric(chr)
  pos <- as.numeric(pos)
  ref <- as.character(ref)
  alt <- as.character(alt)
  N <- as.numeric(N)
  if(length(beta) > 1) beta <- as.numeric(beta)
  if(length(pval) > 1) pval <- as.numeric(pval)
  if(length(af) >1) af <- as.numeric(af)

  #set ids
  id <- ifelse(is.na(id) | id == ".",
               paste(chr[is.na(id) | id == "."], pos[is.na(id) | id == "."], ref[is.na(id) | id == "."], alt[is.na(id) | id == "."], sep=":"), id)

  data <- list(chr=chr, pos=pos, id=id, ref=ref, alt=alt, N=N, af=af, beta=beta, se=se, pval=pval, build=build)
  attr(data, "class") <- "gData"
  return(data)
}

#' is object class gData
#'
#' Checks object for class gData
#' @param x object
#' @export
is.gData <- function(x){
  if(class(x) == "gData"){TRUE}else{FALSE}
}

#' @export
print.gData <- function(x, ...){
  cat(paste(paste0("'",deparse(substitute(x)), "'"), "is an object of class", class(x), "consisiting of", paste(names(x), collapse=",")), "\n")
  cat(paste("data consists of", length(x$chr), "variants", "\n"))
  cat("'chr' is a vector of chromosome ids \n")
  cat("'pos' is a vector of phsyical positions \n")
  cat("'id' is a vector of snp IDs \n")
  cat("'ref' is the vector of reference or A2 alleles \n")
  cat("'alt' is the vector of alternate or A1 alleles \n")
  cat("'N' is the sample size vector \n")
  cat("'AF' is the vector of allele frequencies \n")
  cat("'beta' is a vector of GWAS effect sizes; NA if data source is not GWAS sumstats \n")
  cat("'se' is a vector of standard errors of the beta; NA if data source is not GWAS sumstats \n")
  cat("'pval' is a vector of GWAS P-values; NA if data source is not GWAS sumstats \n")
  cat("'build' is a charater vector of length 1 indicating genome build (hg19 or hg38) \n")
  if(length(x$beta) <2 | length(x$pval) < 2){
    cat("beta and pval vectors are empty; genotype data \n")
  }else{
    cat("GWAS summary statistics data \n")
  }
}

#' Show object gData
#' @param x object of class gData
#' @export
show.gData <- function(x){
  print.gData(x)
}

#' @export
as.data.frame.gData <- function(x, row.names=NULL, optional=FALSE, ..., stringsAsFactors=FALSE){
  if(length(x$beta) >1 & length(x$pval)>1){
    dat <- data.frame(x$chr,x$pos,x$id,x$ref,x$alt, x$beta, x$se, x$pval, x$N, x$af, stringsAsFactors=stringsAsFactors)
    colnames(dat) <- c("CHROM", "POS", "ID", "REF", "ALT", "BETA", "SE", "PVAL", "N", "AF")
  }else{
    dat <- data.frame(x$chr,x$pos,x$id,x$ref,x$alt, x$N, x$af, stringsAsFactors=stringsAsFactors)
    colnames(dat) <- c("CHROM", "POS", "ID", "REF", "ALT", "N", "AF")
  }
  return(dat)
}

#' @export
subset.gData <- function(x, ..., subset){

  x$chr <- x$chr[subset]
  x$pos <- x$pos[subset]
  x$id <- x$id[subset]
  x$ref <- x$ref[subset]
  x$alt <- x$alt[subset]
  if(length(x$beta) > 1) x$beta <- x$beta[subset]
  if(length(x$beta) > 1) x$pval <- x$pval[subset]
  return(x)

}
