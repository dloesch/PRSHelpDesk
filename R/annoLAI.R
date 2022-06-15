#' Annotate VCF with LAI
#'
#' Annotates VCF with local ancestry information by appending results from msp file (LAI) with VCF genotpe records.
#' Can be paired with laiPRS() to estimate a local-ancestry aware PRS
#' @importFrom data.table fread fwrite
#' @param vcf.file path to vcf file. Note: vcf file should just contain phased genotypes
#' @param msp.file path to msp file from rfmix
#' @param chr chromsoome number
#' @param out.prefix prefix for output vcf
#' @param compress compress file? Defaults to true
#' @param bgzip path to utility to compress file. Defaults to "bgzip"
#' @export

annotLAI <- function(vcf.file, msp.file, chr, out.prefix="output", compress=TRUE, bgzip="bgzip"){
  #read in
  vcf <- fread(vcf.file, data.table = FALSE, colClasses = "character")
  msp <- fread(msp.file, data.table=FALSE, colClasses = "character")

  #print out summary
  print(paste("VCF file contains", nrow(vcf), "variants"))
  print(paste("msp file contains", nrow(msp), "ancestry blocks"))

  #prep vcf
  vcf[,2] <- as.numeric(vcf[,2])
  vcf.info <- vcf[1:9]
  vcf <- vcf[10:ncol(vcf)]

  #prep msp
  msp$spos <- as.numeric(msp$spos)
  msp$epos <- as.numeric(msp$epos)

  #subject ids
  ids <- colnames(vcf)
  hap1 <- paste0(ids, ".0")
  hap2 <- paste0(ids, ".1")

  #collapse haplotpyes to per-sample. If msp file indicates a 1 and a 2, results in 1|2
  foo <- msp[1:6]
  for(id in ids){
    hap1 <- paste0(id, ".0")
    hap2 <- paste0(id, ".1")
    foo[[id]] <- paste0(msp[[hap1]], "|", msp[[hap2]])
  }

  #prepare output
  out.file <- paste0(out.prefix, ".LA.chr", chr, "vcf")
  cat('##fileformat=VCFv4.2\n', file=out.file)
  cat('##FILTER=<ID=PASS,Description="All filters passed">\n', file=out.file, append = TRUE )
  cat('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' , file=out.file, append = TRUE )
  cat('##FORMAT=<ID=LA,Number=1,Type=String,Description="Local Ancestry status">\n' , file=out.file, append = TRUE )
  cat(paste0("##contig=<ID=chr", chr, '>\n') , file=out.file, append = TRUE )
  cat('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">\n' , file=out.file, append = TRUE )
  cat('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">\n' , file=out.file, append = TRUE )

  #extract just LAI matrix
  temp <- as.matrix(foo[7:ncol(foo)])
  if(ncol(temp) != ncol(vcf) ) stop("vcf and lai sample size do not match")

  for(i in 1:nrow(foo)){
    if(i == nrow(foo)){
      temp.vcf <- as.matrix(vcf[vcf.info$POS >= foo$spos[i] & vcf.info$POS <= foo$epos[i],])
      temp.info <- vcf.info[vcf.info$POS >= foo$spos[i] & vcf.info$POS <= foo$epos[i],]
    }else{
      temp.vcf <- as.matrix(vcf[vcf.info$POS >= foo$spos[i] & vcf.info$POS < foo$epos[i],])
      temp.info <- vcf.info[vcf.info$POS >= foo$spos[i] & vcf.info$POS < foo$epos[i],]
    }

    if(nrow(temp.vcf) == 0) next
    #make msp file same dimensions as temp.vcf
    temp.msp <- matrix(temp[i,], nrow = 1)
    for(j in 1:(nrow(temp.vcf)-1)) temp.msp <- rbind(temp.msp, matrix(temp[i,], nrow=1))

    vcf.chunk <- paste(temp.vcf,temp.msp, sep=":")
    vcf.chunk <- matrix(data = vcf.chunk, nrow = nrow(temp.vcf), ncol = ncol(temp.vcf))

    #prep vcf file
    new.vcf <- as.data.frame(vcf.chunk, stringsAsFactors = FALSE)
    colnames(new.vcf) <- colnames(temp)
    new.vcf <- cbind(temp.info, new.vcf)
    new.vcf$FORMAT <- paste("GT", "LA", sep=":")

    #write-out vcf
    if(i == 1){
      fwrite(new.vcf, out.file, append = TRUE, sep='\t', col.names = TRUE, row.names = FALSE, quote=FALSE)
    }else{
      fwrite(new.vcf, out.file, append = TRUE, sep='\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
    }
  }

  if(compress == TRUE){
    system2(bgzip, out.file)
  }
  print("done")
}
