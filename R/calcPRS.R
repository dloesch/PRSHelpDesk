#' Calculate polygenic risk score
#'
#' Calculates polygenic risk score with snps present in mgData object (merged target genotype + GWAS sumstats)
#' Can be used to validate pre-trained PRS. Also implements a simple pruning/thresholding procedure using a windowing procedure.
#' @param mgData GWAS summary statistics, can be data.frame or filepath
#' @param vcfData File path to VCF, BCF, BED
#' @param imputed Is data imputed? Defaults to FALSE
#' @param CGAT Keep CG/AT sites? Defaults to FALSE
#' @param AF filter by allele frequencies? Defaults to FALSE
#' @param af.cutoff Allele frequency difference cutoff. Defaults to NA
#' @param prune Prune by kb window? defaults to false
#' @param kb kb window for prune
#' @param thresh filter by pvalue threshold? Defaults to FALSE
#' @param pval pvalue threshold. Defaults to 5E-08
#' @return returns data.frame with columns ID, PRS_RAW, PRS_AVG, PRS_SCALED. If save is TRUE, also outputs
#' text file with model SNPs and a VCF file of SNPs included in model.
#' @export
#' @importFrom data.table fread fwrite
#' @examples
#'  \dontrun{
#'  calcPRS(sumstats, "data.vcf.gz", prefix="cohort", save=TRUE)
#'  }
calcPRS <- function(mgData, vcfData, imputed=FALSE, CGAT=FALSE, AF=FALSE, af.cutoff=NA, prune=FALSE, kb=1000,
                    thresh=FALSE, pval=5E-08){

  x <- mgData

  if(CGAT == FALSE){
    x <- subset(x, x$cgat == FALSE)
  }

  #filter
  x <- subset(x, x$filter == "PASS")

  #filter by af, if desired
  if(AF == TRUE){
    af.filter <- testAF(x$af1, x$af2, cutoff = af.cutoff, align = x$info)
    x <- subset(x, af.filter == "PASS")
  }

  #prune, if desired
  if(prune == TRUE){
    index.snps <- c()
    for(chr in unique(x$chr)){
      start <- min(x$pos[x$chr == chr])
      end <- max(x$pos[x$chr == chr])
      i <- start
      j <- start + kb*1000
      while(i <= end){
        min.pval <- min(x$pval[x$pos >= i & x$pos <= j  & x$chr == chr])
        index.snps <- c(index.snps, x$id1[x$pval == min.pval & x$chr == chr & x$pos >= i & x$pos <= j][1])
        i <- i + (kb*1000)*0.9
        j <- ifelse(j > end, end, i + kb*1000)
      }

    }
    #remove duplicates that might occur due to the overlapping windows
    index.snps <- unique(index.snps)
    x <- subset(x, x$id1 %in% index.snps)
  }

  #threshold by pvalue, if desired
  if(thresh == TRUE){
    x <- subset(x, x$pval <= pval)
  }

  #subset vcf file

  vcf <- subset(vcfData, rownames(vcfData$geno) %in% x$id1)

  #prepare betas
  betas <- ifelse(x$info == "SWITCH" | x$info == "FLIP_SWITCH", x$beta*-1, x$beta)

  #prepare genotype data
  geno <- vcf$geno

  if(imputed == TRUE){
    geno <- imputed_dosage(geno)
  }else{
    geno <- to_dosage(geno)
  }
  geno <- geno*betas
  prs <- unlist(unname(colSums(geno, na.rm = TRUE)))

  data <- data.frame(ID=colnames(geno),
                     NSNPS=length(betas),
                     PRS_RAW=prs,
                     PRS_AVG=prs/(length(betas)*2),
                     PRS_SCALE=scale(prs),
                     stringsAsFactors = FALSE)

  return(data)
}
