#' Clump+Thresh PRS
#'
#' Uses PLINK to obtain a P+T/C+T set of summary sttatitics given an LD reference file, target genotype file, GWAS sumstats, R2, kb window, and pvalue
#' parameters. These paramters can be learned by writing a loop or using PRSice-2, or a combination of the two. The resulting pruned
#' file can then be used with calcPRS2() to estimate a PRS.
#' @importFrom data.table fwrite fread
#' @importFrom stats complete.cases
#' @param target.file prefix of target plink fileset
#' @param ref.file prefix of LD ref plink fileset
#' @param sumstats.file path to summary stats file
#' @param snp.col name of SNP column
#' @param chr.col chromosome column
#' @param bp.col pos/basepair column
#' @param a1.col A1/effect allele/alt allele column
#' @param beta.col beta column
#' @param pval.col pvalue column
#' @param plink path to plink
#' @param kb size of LD window. Default is 250 kb
#' @param R2 R2 parameter for LD. Default is 0.1
#' @param pval P-value parameter for thresholding. Default is 5E-08
#' @export
#'

clumpPRS <- function(target.file, ref.file, sumstats.file, snp.col, chr.col, bp.col, a1.col, beta.col, pval.col, plink="plink",
                     kb=250, R2=0.1, pval=5E-08){

  target.bim <- fread(paste0(target.file, ".bim"), data.table = FALSE)

  snps <- target.bim[,2]

  sumstats <- fread(sumstats.file, data.table = FALSE)

  sumstats <- sumstats[c(snp.col, chr.col, bp.col, a1.col, beta.col, pval.col)]
  colnames(sumstats) <- c("SNP","CHR", "BP", "A1", "BETA", "PVAL")
  sumstats <- sumstats[complete.cases(sumstats),]
  #only keep snps present in target dataset
  sumstats <- sumstats[sumstats$SNP %in% snps,]

  #make it in plink format and write out
  #CHR	SNP BP A1 TEST NMISS BETA/OR STAT P
  plink.format <- data.frame(CHR=sumstats$CHR,
                      SNP=sumstats$SNP,
                      BP=sumstats$BP,
                      A1=sumstats$A1,
                      NMISS=rep(0, times=nrow(sumstats)),
                      BETA=sumstats$BETA,
                      STAT=rep(0, times=nrow(sumstats)),
                      P=sumstats$PVAL,
                      stringsAsFactors = FALSE)

  out.file <- tempfile()
  fwrite(plink.format, paste0(out.file, ".txt"), col.names = TRUE, row.names = FALSE, sep='\t', quote=FALSE, na = "NA")
  system2(plink,
          paste("--bfile", ref.file, "--clump", paste0(out.file, ".txt"), "--clump-kb", kb, "--clump-r2", R2, "--clump-p1", pval, "--clump-p2", pval,
                "--out", out.file))

  clump <- fread(paste0(out.file, ".clumped"), data.table = FALSE)
  sumstats$BETA <- ifelse(sumstats$SNP %in% clump$SNP, sumstats$BETA, 0)

  #filter by pval threshold
  sumstats <- sumstats[sumstats$PVAL <= pval,]
  sumstats <- sumstats[c("SNP", "A1", "BETA")]
  return(sumstats)

}
