#' File type of genetic data
#'
#' Infers filetype from file extension
#' @importFrom tools file_ext
#' @param x filename
#' @export
file.type <- function(x){
  file.type <- ifelse(file_ext(x)=="bcf", "bcf",
                      ifelse(file_ext(x) == "bed", "bed",
                             ifelse(file_ext(sub(".gz", "", x)) == "vcf", "vcf", "other")))
  if(file.type == "other"){
    ifelse(file.exists(paste0(x, ".bed")), "bed", "other")
  }
  return(file.type)
}
