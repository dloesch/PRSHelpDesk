#' Concatenate data split by chromosome
#'
#' Concatenates data split by chromosome, i.e. GWAS data.
#' Provide fie name with a space for chromosome number.
#' @param filename base name of file
#' @return data.frame with concatenated data
#' @export
#' @importFrom data.table fread
#' @examples
#' \dontrun{
#' concatData("GWAS.results.chr .txt")
#' }
concatData <- function(filename){

  part1 <- unlist(unname(strsplit(filename, split=" ")))[1]
  part2 <- unlist(unname(strsplit(filename, split=" ")))[2]
  for(chr in 1:22){
    data.file <- paste0(part1, chr, part2)
    foo <- fread(data.file, data.table=FALSE)
    if(chr ==1){
      data <- foo
    }else{
      data <- rbind(data, foo)
    }
  }
  return(data)
}
