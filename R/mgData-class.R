#' merge objects of class gData
#'
#' S3 method for merging class gData and contructor function for class mgData
#' @aliases mgData
#' @param x gData object to merge
#' @param y gData object to merge
#' @param ... further arguments passed to or from other methods.
#' @param by column for merging. Default is by chrom and pos
#' @param all.x Defaults to true
#' @param sort Defaults to FALSE
#' @export
#' @return object of class mgData
#'
merge.gData <- function(x,y,...,by=c("CHROM", "POS"), all.x=TRUE, sort=FALSE){
  target <- x
  gwas <- y
  if(target$build != gwas$build) stop("genome builds do not match")
  build <- target$build

  #parse first gData from target
  x <- as.data.frame(target)
  colnames(x)[3:5] <- paste0(colnames(x)[3:5], "1")

  #parse second gData file from gwas
  y <- as.data.frame(gwas)
  colnames(y)[3:5] <- paste0(colnames(y)[3:5], "2")
  z <- merge(x, y, by=by, all.x=all.x,  sort=sort)

  #identify CGAT
  cgat <- mapply(FUN = findCGAT, z$REF1, z$ALT1, SIMPLIFY = TRUE)
  #match alleles
  info <- mapply(FUN= check_alleles, z$REF1, z$ALT1, z$REF2, z$ALT2, SIMPLIFY = TRUE)

  #filter out allele mismatches and alleles not in reference
  filter <- ifelse(info == "MISMATCH" | info == "NOT_IN_REFERENCE", "FAIL", "PASS")
  #flag duplicates
  dupes <- duplicated(z$ID2)
  info <- ifelse(dupes, "DUPLICATE", info)
  filter <- ifelse(dupes, "FAIL", filter)

  m.gData <- list(chr=z$CHROM, pos=as.numeric(z$POS), id1=as.character(z$ID1), ref1=as.character(z$REF1), alt1=as.character(z$ALT1), af1=z$AF.x, n1=z$N.x,
                  id2=as.character(z$ID2), ref2=as.character(z$REF2), alt2=as.character(z$ALT2), af2=z$AF.y, n2=z$N.y, beta=z$BETA, pval=z$PVAL,
                  build=build, cgat=cgat, info=info, filter=filter)

  attr(m.gData, "class") <- "mgData"
  return(m.gData)
}

#' @export
as.data.frame.mgData <- function(x, row.names=NULL, optional=FALSE, ..., stringsAsFactors = FALSE){
  dat <- data.frame(x$chr, x$pos, x$id1, x$ref1, x$alt1, x$af1, x$n1,
                      x$id2, x$ref2, x$alt2, x$af2, x$n2,
                      x$cgat, x$info, x$filter, x$beta, x$pval, stringsAsFactors = stringsAsFactors)
  colnames(dat) <- c("#CHROM","POS" ,"ID1" , "REF1", "ALT1" , "AF1", "N1",
                       "ID2", "REF2", "ALT2", "AF2", "N2",
                       "CGAT", "INFO", "FILTER", "BETA", "PVAL")

  return(dat)
}

#' Check of class mgData
#' Simple check if class mgData
#' @param x object
#' @export
is.mgData <- function(x){
  if(class(x) == "mgData"){TRUE}else{FALSE}
}

#' @export
print.mgData <- function(x, ...){
  cat(paste(paste0("'",deparse(substitute(x)), "'"), "is an object of class", class(x), "consisiting of", paste(names(x), collapse=",")), "\n")
  cat(paste("data consists of", length(x$chr), "variants after merging two gData objects", "\n"))
  cat("'chr' is a vector of chromosome ids \n")
  cat("'pos' is a vector of phsyical positions \n")
  cat("'id1' and 'id2' are vectors of snp IDs \n")
  cat("'ref1' and 'ref2' are vectors of reference or A2 alleles \n")
  cat("'alt1' and 'alt2; are vectors of alternate or A1 alleles \n")
  cat("'beta' is a vector of GWAS effect sizesl NA if netiher data source was GWAS sumstats \n")
  cat("'pval' is a vector of GWAS P-values; NA if netiher data source was GWAS sumstats \n")
  cat("'build' is a charater vector of length 1 indicating genome build (hg19 or hg38) \n")
}

#' show object of class mgData
#' @param x object of class mgData
#' @export
show.mgData <- function(x){
  print.mgData(x)
}

#' @export
subset.mgData <- function(x, subset, ...){

  x$chr <- x$chr[subset]
  x$pos <- x$pos[subset]
  x$id1 <- x$id1[subset]
  x$ref1 <- x$ref1[subset]
  x$alt1 <- x$alt1[subset]
  x$af1 <- x$af1[subset]
  x$n1 <- x$n1[subset]
  x$id2 <- x$id2[subset]
  x$ref2 <- x$ref2[subset]
  x$alt2 <- x$alt2[subset]
  x$af2 <- x$af2[subset]
  x$n2 <- x$n2[subset]
  x$cgat <- x$cgat[subset]
  x$info <- x$info[subset]
  x$filter <- x$filter[subset]
  x$beta <- x$beta[subset]
  x$pval <- x$pval[subset]
  return(x)

}
