#' find CG/AT sites
#'
#' Flags CG/AT sites as they are strand ambiguous
#' @param ref ref allele
#' @param alt alt allele
#' @return TRUE or FALSE
#' @export
#' @examples
#' findCGAT("C", "G")
findCGAT <- function(ref,alt){
  ref <- toupper(ref)
  alt <- toupper(alt)
  CGAT_filter <- ifelse(ref == "C" & alt == "G", TRUE, ifelse(ref == "G" & alt == "C", TRUE,
                                                              ifelse(ref == "A" & alt == "T", TRUE,
                                                                     ifelse(ref == "T" & alt == "A", TRUE, FALSE))))
  return(CGAT_filter)
}


