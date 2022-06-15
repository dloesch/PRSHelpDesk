#' convert genotype files
#'
#' Wrapper function for handling file conversions
#' @param geno.file filepath to genotype data
#' @param from Filetype of original genotype data
#' @param to Desired filetype
#' @param prefix Desired prefix. Defaults to output
#' @export
convertGENO <- function(geno.file,from=c("bed", "bcf", "vcf"),to=c("bed", "bcf", "vcf"), prefix="output"){
  if(!file.exists(geno.file)){
    stop("error: file does not exist")
  }

  #check to/from arguments:###
  from <- tolower(from)
  to <- tolower(from)
  if(length(to) > 1){
    stop("error: select one file type")
  }
  if(length(from) > 1){
    stop("error: select one file type")
  }
  if(from == to){stop("error: to and from are same file type")}
  if(!from %in% c("bed", "bcf", "vcf")){ stop("error: file type needs to be bed, bcf, vcf")}
  if(!to %in% c("bed", "bcf", "vcf")){ stop("error: file type needs to be bed, bcf, vcf")}

  ##done##

  if(from == "bed"){
    geno.file <- sub(".bed", "", geno.file)
    print("converting to vcf")
    system2("plink", paste("--bfile", geno.file, "--recode vcf-iid bgz --out", prefix), stdout = FALSE)
  }else{

    if(to == "bed"){
      if(from == "bcf"){
        system2("plink", paste("--bcf", geno.file, "--make-bed --out", prefix), stdout = FALSE)
      }else{
        system2("plink", paste("--vcf", geno.file, "--make-bed --out", prefix), stdout = FALSE)
      }
    }else{
      if(from== "bcf"){
        system2("bcftools", paste("view -Oz -o", paste0(prefix, ".vcf.gz"), geno.file))
      }else{
        system2("bcftools", paste("view -Ob -o", paste0(prefix, ".bcf"), geno.file))
      }
    }
  }

}
