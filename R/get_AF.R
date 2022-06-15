#' get allele frequencies
#'
#' calculates allele frequencies external to R using appropriate software.
#' For VCF/BCF files, will first check to see if INFO contains AF
#' @param geno.file file path to genotype file
#' @param plink path to plink. Default assumes availavle in path
#' @param bcftools path to bcftools. Default assumes available in path
#' @return data.frame with columns SNP and AF
#' @importFrom utils read.table
#' @importFrom data.table fread
#' @export

get_AF <- function(geno.file, plink="plink", bcftools="bcftools"){
  if(!file.exists(geno.file)) stop("error: unable to read file")

  #check to see if allele frequencies are already calculated in vcf/bcf file
  af <- checkAF(geno.file)
  if(af){
    p <- pipe(paste(bcftools,  "query -f '%ID %INFO/AF %INFO/AN\n' ", geno.file), "rt")
    af <- read.table(p, stringsAsFactors = FALSE)
    close(p)
    colnames(af) <- c("SNP", "AF", "AN")
  }else{
    #use appropriate tool to calculate allele frequencies
    if(file.type(geno.file) == "bed"){
      input <- paste("--bfile", sub(".bed", "", geno.file))
    }
    else if(file.type(geno.file) == "vcf"){
      input <- paste("--vcf", geno.file)
    }
    else if(file.type(geno.file) == "bcf"){
      input <- paste("--bcf", geno.file)
    }else{
      stop("Supports bcf/vcf/plink files")
    }
      temp <- tempfile()
      system2(plink, paste(input, "--freq --keep-allele-order  --out", temp), stdout = FALSE)
      af <- fread(paste0(temp, ".frq"), data.table=FALSE)
      af <- af[c("SNP", "MAF", "NCHROBS")]
      colnames(af)[2:3] <- c("AF", "AN")
      system2("rm", paste0(temp,"*"))
  }
  return(af)
}
