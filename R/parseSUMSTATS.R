#' Parse GWAS sumstats
#'
#' Parse GWAS sumstats and return as gData object
#' @importFrom data.table fread fwrite
#' @importFrom stats complete.cases
#' @param sumstats.file Sumstats file
#' @param snp.col column name for snp IDs.
#' @param beta.col column name for betas.
#' @param pval.col column name for pvals.
#' @param chrom.col column name for chromosome number.
#' @param pos.col column name for physical position.
#' @param ref.col column name for reference/other allele.
#' @param alt.col column name for alt/effect allele.
#' @param af.col column name for effect allele frequency. Defaults to NA
#' @param se.col column name for standard error of the beta. Defaults to NA
#' @param N.col column name for sample size. defaults to NA
#' @param N sample size, if single value. Defaults to NA
#' @param OR effect sizes are odds ratios; convert to betas. Defaults to false
#' @param build Genome build. defaults to hg19
#' @param bim Bim file containing rsids, if needed. Defaults to NA
#' @export

parseSUMSTATS <- function(sumstats.file,
                      snp.col, beta.col, pval.col, chrom.col, pos.col, ref.col, alt.col,
                      af.col=NA, se.col=NA, N.col=NA, N=NA,
                      OR=FALSE, build="hg19", bim=NA){

  #check if sumstats file exists
  if(!file.exists(sumstats.file)) stop("error: unable to read GWAS sumstats")

  data <- fread(sumstats.file, data.table = FALSE)

  ###prep input###
  #only keep "chr" prefix if build38
  data[[chrom.col]] <- gsub("chr", "", data[[chrom.col]])
  data[[beta.col]] <- as.numeric(data[[beta.col]])
  data[[pval.col]] <- as.numeric(data[[pval.col]])
  data[[pos.col]] <- as.numeric(data[[pos.col]])
  data[[ref.col]] <- toupper(data[[ref.col]])
  data[[alt.col]] <- toupper(data[[alt.col]])
  if(is.na(af.col)) {
    af.col <- "AF"
    data$AF <- NA
  }
  if(is.na(se.col)){
    se.col <- "SE"
    data$SE <- NA
  }
  #specify study size
  if(is.na(N.col)){
    N.col <- "N"
    data$N <- NA
  }
  if(!is.na(N)) data$N <- N

  #convert odds ratios to betas
  if(OR == TRUE){
    data[[beta.col]] <- log(data[[beta.col]])
  }

  #make temporary snp id column, if necessary
  if(snp.col == "NONE" | snp.col == "NA" | is.na(snp.col)){
    data$SNP <- paste(data[[chrom.col]], data[[pos.col]], data[[ref.col]], data[[alt.col]], sep=":")
    snp.col <- "SNP"
  }


  #organize data
  data <- data[c(snp.col, chrom.col, pos.col, ref.col, alt.col, beta.col, pval.col, se.col, af.col, N.col)]
  colnames(data) <- c("SNP", "CHR", "BP", "A2", "A1", "BETA", "P", "SE", "AF", "N")

  #sort
  data <- data[order(data$CHR, data$BP),]
  #add "chr" prefix back for build 38
  if(build == "hg38" | build == "HG38" | build == "GC38"){
    data$CHR <- paste0("chr", data$CHR)
  }

  #exclude NA sites
  data <- data[complete.cases(data$SNP),]

  #make sure duplicated snps are not in the data
  #exclude them all, not just keep first occurance
  dupes <- data$SNP[duplicated(data$SNP)]
  data <- data[!data$SNP %in% dupes,]

  #print warning if GWAS data does not contain rsids
  if(length(grep("rs", data$SNP)) < 0.8*nrow(data)){
    print("Warning: less than 80% of the data contains rsids. Provide bim file with rsids to fix")
  }


  if(!is.na(bim)){
    if(!file.exists(bim)) stop("error: unable to read bim file")

    bim <- fread(bim, data.table = FALSE)
    colnames(bim) <- c("CHR", "SNP.2", "CM", "BP", "A1.2", "A2.2")
    data <- merge(data, bim, by=c("CHR", "BP"), all=FALSE, sort=FALSE)

    data$CHECK <- mapply(check_alleles, data$A2, data$A1, data$A2.2, data$a1.2)
    data <- data[data$CHECK !=  "ALLELE_MISMATCH" & data$CHECK != "NOT_IN_REFERENCE",]
    data$SNP <- data$SNP.2
  }

  #prep to return data
  sumstats <- data[c("SNP", "CHR", "BP", "A2", "A1", "BETA", "P", "SE", "AF", "N")]

  data <- gData(chr = sumstats$CHR, pos = sumstats$BP, id = sumstats$SNP, ref = sumstats$A2, alt = sumstats$A1, N = sumstats$N,
                beta = sumstats$BETA, se = sumstats$SE, pval = sumstats$P, af = sumstats$AF, build = build)


  return(data)

}
