#' check alleles
#'
#' Checks if alleles match, are strand flipped, switched (ref/alt swap)
#' @param ref ref allele for target data
#' @param alt allele for target data
#' @param r.ref ref allele for reference data
#' @param r.alt alt allele for reference data
#' @return genotpe status based on allele labels
#' @export
#' @examples
#' check_alleles("A", "T", "T", "A")
check_alleles <- function(ref, alt, r.ref, r.alt){
  ref <- toupper(ref)
  alt <- toupper(alt)
  r.ref <- toupper(r.ref)
  r.alt <- toupper(r.alt)

  labels <- c("NOT_IN_REFERENCE", "SAME", "SWITCH", "FLIP", "FLIP_SWITCH", "MISMATCH")

  if(is.na(r.ref) | is.na(r.alt)){
    label <- labels[1]
  }else{
    label <- ifelse(r.ref == ref & r.alt == alt, labels[2],
                    ifelse(r.ref == alt & r.alt == ref, labels[3],
                           ifelse(r.ref == flip_allele(ref) & r.alt == flip_allele(alt), labels[4],
                                  ifelse(r.ref == flip_allele(alt) & r.alt == flip_allele(ref), labels[5], labels[6]))))
  }
  return(label)

}
