#' class vcfData
#'
#' Stores just the genotype data as a data.frame and calculates frequencies if AF=TRUE
#' @param vcf data.frame of full data when reading vcf into R
#' @export

vcfData <- function(vcf){
  #run some checks
  if(!is.data.frame(vcf)){
    tryCatch({
      vcf <- as.data.frame(vcf)
    },
    error=function(e)
    print("error:input needs to be data.frame in vcf format")
    )
  }
  #needs to have at least one genotype record in one site, i.e. 10 columns and 1 row
  if(ncol(vcf) < 10 | nrow(vcf)<1) stop("error: input needs to be data.frame in vcf format")

  colnames(vcf)[1:5] <- c("CHROM", "POS", "ID", "REF", "ALT")
  geno <- as.data.frame(vcf[10:ncol(vcf)])
  N <- ncol(geno)
  M <- nrow(geno)

  #set ids
  id <- ifelse(is.na(vcf$ID) | vcf$ID == ".", paste(vcf$`#CHROM`, vcf$POS, vcf$REF, vcf$ALT, sep=":"), vcf$ID)
  rownames(geno) <- id

  data <- list(geno=geno, N=N, M=M)
  attr(data, "class") <- "vcfData"
  return(data)
}

#' @export
as.data.frame.vcfData <- function(x, row.names=NULL, optional=FALSE, ...){
  snp.ids <- rownames(x$geno)
  ids <- colnames(x$geno)
  rownames(x$geno) <- NULL
  dat <- data.frame(snp.ids, x$geno)
  colnames(dat) <- c("ID", ids)

  return(dat)
}

#' Simple check if object is vcfData
#' @param x object
#' @export
is.vcfData <- function(x){
  if(class(x) == "vcfData"){TRUE}else{FALSE}
}

#' @export
print.vcfData <- function(x, ...){
  cat(paste(paste0("'",deparse(substitute(x)), "'"), "is an object of class", class(x), "consisiting of", paste(names(x), collapse=",")), "\n")
  cat(paste("'geno' is a data.frame of", x$M, "markers (M) and", x$N, "subjects (N)", "\n"))
}

#' show object of class vcfData
#' @param x object of class vcfData
#' @export
show.vcfData <- function(x){
  print(x)
}

#' @export
subset.vcfData <- function(x,subset, by="marker", ...){
  if(by == "marker"){
    x$geno <- x$geno[subset,]
    x$M <- nrow(x$geno)
  }else{
    x$geno <- x$geno[,subset]
    x$N <- ncol(x$geno)
  }
  return(x)
}

#' Concatenate vcfData objects
#'
#' concatenate vcfData objects. Checks if number of subjects is the same. If both
#' have allele frequencies, they will be combined via c()
#'
#' @param x vcfData object 1
#' @param y vcfData object 2
#' @return vcfData object
#' @export
concat.vcfData <- function(x, y){
  if(x$N == y$N){
    x$geno <- rbind(x$geno, y$geno)
    x$M <- x$M + y$M
    return(x)
  }else{
    stop("error: different number of subjects ")
  }

}
