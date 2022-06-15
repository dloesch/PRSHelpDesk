#' Calculate PRS, version 2
#'
#' Estimates a PRS using PLINK --score function given paths to genotype file and sumstats file,
#' i.e. file of posterior effect sizes from PRS-CS
#'
#' @importFrom data.table fwrite fread
#' @importFrom stats complete.cases
#' @param geno.file path to genotype file
#' @param sumstats.file path to sumstats file
#' @param snp.col name of snp column
#' @param a1.col name of a1 column
#' @param beta.col name of beta column
#' @param plink path to plink
#' @export

calcPRS2 <- function(geno.file, sumstats.file, snp.col, a1.col, beta.col, plink="plink"){

  sumstats <- fread(sumstats.file, data.table = FALSE, na.strings = "")

  sumstats <- sumstats[c(snp.col, a1.col, beta.col)]
  sumstats <- sumstats[complete.cases(sumstats),]

  temp.stats <- tempfile()
  fwrite(sumstats, temp.stats, col.names = FALSE, row.names = FALSE, sep='\t', quote=FALSE)

  temp.prs <- tempfile()

  if(file.type(geno.file) == "vcf"){
    input <- paste("--vcf", geno.file)
  }else if( file.type(geno.file) == "bcf"){
    input <- paste("--vcf", geno.file)
  }else{
    input <- paste("--bfile", sub(".bed", geno.file))
  }

  system2(plink,
          paste(input, "--score", temp.stats, "1 2 3 sum --out", temp.prs))
  foo <- fread(paste0(temp.prs, ".profile"), data.table = FALSE)
  prs <- data.frame(FID=foo$FID, IID=foo$IID, NSNPS=foo$CNT/2, PRS_RAW=foo$SCORESUM,
                    stringsAsFactors = FALSE)

  prs$PRS_AVG <- prs$PRS_RAW/(prs$NSNPS*2)
  prs$PRS_SCALE <- scale(prs$PRS_RAW)

  return(prs)
}
