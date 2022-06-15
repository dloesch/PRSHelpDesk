#' parse genotype files
#'
#' parses genotype files, converting either to class gData if by-site information is wanted
#' or class vcfData if genotype data is wanted
#' @param geno.file file path to genotype data
#' @param build genome build. Either hg19 or hg38. Defaults to hg19
#' @param geno Read in genotype data? Defaults to FALSE
#' @param both Read in full data? Defaults to FALSE
#' @param af Calculate allele frequencies? Deafuls to FALSE
#' @param plink path to plink
#' @param bcftools path to bcftools
#' @export
#' @importFrom utils read.table
#' @importFrom data.table fread
parseVCF <- function(geno.file, build="hg19", geno=FALSE, both=FALSE, af=FALSE, plink="plink", bcftools="bcftools"){
  #check inputs
  if(!file.exists(geno.file)) stop("error: unable to read file")
  build <- tolower(build)
  if(build != "hg19" & build != "hg38") stop("error: build needs to be hg19 or hg38")

  #extract genotype if desired
  if(geno == TRUE & both == FALSE){
    vcf <- read_vcf(geno.file)
    data <- vcfData(vcf)
    return(data)
  }else if(both == TRUE){
    vcf <- read_vcf(geno.file)
    N <- ncol(vcf) -9
    if(af == TRUE){
      af_dat <- get_AF(geno.file, plink = plink, bcftools = bcftools)
      af_dat <- af_dat$AF
    }else{
      af_dat <- NA
    }
    gdata <- gData(chr = vcf[,1], pos = vcf[,2], id = vcf[,3], ref = vcf[,4], alt = vcf[,5], N=N, af = af_dat, build = build)
    vcf <- vcfData(vcf)
    data <- list(gdata=gdata, vcf=vcf)
    return(data)
  }else{
    #does not read in genotype information
    if(file.type(geno.file) == "bed"){
      geno.file <- sub(".bed", ".bim", geno.file)
      data <- fread(geno.file, data.table=FALSE)
      data <- data[c(1,4,2,6,5)]
    }else if(file.type(geno.file) != "other"){
      #pipe in. will test if writing to temp file and mult-threaded reading is faster
      p <- pipe(paste(bcftools, "query -f '%CHROM  %POS %ID %REF  %ALT{0}\n'", geno.file), "rt")
      data <- read.table(p, stringsAsFactors = FALSE)
      close(p)
    }else{
      stop("file needs to be BCF/VCF/BED")
    }
    if(af == TRUE){
      af_dat <- get_AF(geno.file, plink = plink, bcftools = bcftools)
      N <- round(mean(af_dat$AN/2))
      af_dat <- af_dat$AF
    }else{
      af_dat <- NA
      N <- NA
    }
    gdata <- gData(chr = data[,1], pos = data[,2], id = data[,3], ref = data[,4], alt = data[,5], N=N, af= af_dat, build = build)
  }

}
