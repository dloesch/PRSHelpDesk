#'compress file
#'
#'compresses file using system call to gzip or bgzip.
#'bgzip needs to be installed and in path. Otherwise, use gzip
#' @param file.name Path to file for compression
#' @param bgzip Use bgzip.
#' @return zipped file
#' @export
#' @examples
#' \dontrun{
#' zip_file("data.vcf", bgzip=TRUE)
#' }
zip_file <- function(file.name, bgzip=TRUE){
  if(bgzip == TRUE){
    system2("which", "bgzip") -> test
    if(test == 0){
      print("compressing with bgzip")
      system2("bgzip", paste("-f", file.name))
    }else{
      print("unable to compress with bgzip")
    }
  }else{
    system2("which", "gzip") -> test
    if(test ==0){
      print("compressing with gzip")
      system2("gzip", paste("-f", file.name))
    }else{
      print("unable to compress with gzip")
    }
  }
}
