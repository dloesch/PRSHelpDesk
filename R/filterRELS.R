#' Remove close relatives
#'
#' Simple Algorithm for removing relatives. Provide kinship matrix and algorithm will quickly identify an optimal set of relatives to remove.
#' First, highly connected subjects are determined by counting the number of pairwise comparisons they appear in (after filtering)
#' Each highly connected subject is removed, in order of number of connections, until all pairwise comparisons contain unique subjects
#' The remaining pairwise comparisons are resolved via genotype missingness rate.
#' NOTE: input can be kinship matrix where rownames are subject ids OR three column data.frame of pairwise relationships wich columns
#' ID1, ID2, kinship coefficient
#' @param kin.mat kinship matrix. Can be matrix or pairwise relationships
#' @param pairwise Pairwise or matrix? Defaults to FALSE
#' @param degree 1st,2nd, or 3rd degree threshold. Deafults to 3rd degree.
#' @param plink.file Path to plink fileset
#' @param plink path to plink binary
#' @param FID include FIDs in output? Defaults to TRUE
#' @export
filterRELS <- function(kin.mat, pairwise=FALSE,  degree="3rd", plink.file, plink="plink", FID=TRUE){
  mat <- kin.mat
  thresh <- ifelse(degree == "1st", 0.177, ifelse(degree == "2nd", 0.0884, 0.0442))
  if(pairwise == FALSE){
    #extracts pairwise relationships from matrix
    #filter by desired threshold
    mat[mat < thresh] <- 0

    #remove upper triangle of matrix
    mat[upper.tri(mat)] <- 0

    ids <- rownames(mat)
    for(i in 1:length(ids)){
      id <- ids[i]
      #get vecfor of all comparisons with individual
      temp <- mat[i,]
      #flag relatives (matrix is sparse)
      rels <- ids[temp > 0]
      rels <- rels[rels != id] #remove self
      #convert to data frame
      if(i == 1){
        rel_pairs <- data.frame(ID1=rep(id, times=length(rels)), ID2=rels, stringsAsFactors = FALSE)
      }else{
        foo <- data.frame(ID1=rep(id, times=length(rels)), ID2=rels, stringsAsFactors = FALSE)
        rel_pairs <- rbind(rel_pairs, foo)
      }
    }
  }else{
    #assumes 3 columns: ID1, ID2, kinship coefficient
    #filter by threshold
    mat <- kin.mat[kin.mat[,3] >=thresh,]

    #remove self, if present
    mat <- mat[mat[,1] != mat[,2],]
    #
    rel_pairs <- mat[1:2]
  }

  #resolve relative pairs
  n_times <- c()
  ids <- unique(c(rel_pairs$ID1, rel_pairs$ID2))
  for(i in ids){
    n <- nrow(rel_pairs[rel_pairs$ID1 == i,]) + nrow(rel_pairs[rel_pairs$ID2 == i,])
    n_times <- c(n_times, n)
  }

  #flag subjets who appear multiple times
  counts <- data.frame(ID=ids, N=n_times, stringsAsFactors = FALSE)
  counts <- counts[counts$N > 1,]

  #start filtering, starting with subject who appear multiple times
  temp <- rel_pairs
  temp_drop <- counts$ID
  drop <- c()
  for(i in counts$ID[order(counts$N, decreasing = TRUE)]){

    if(i %in% unique(c(temp$ID1, temp$ID2))){
      drop <- c(drop, i)
      temp <- temp[temp$ID1 != i,]
      temp <- temp[temp$ID2 != i,]
    }else{
      next
    }

  }

  #resolve remaining relative pairs using genotype missingess
  print("calculating genotype missingess rates")
  temp.miss <- tempfile()
  system2(plink, paste("--bfile", plink.file, "--missing --out", temp.miss))

  miss <- data.table::fread(paste0(temp.miss, ".imiss"), data.table = FALSE)
  print(paste("deleting files at", temp.miss))
  system2("rm", paste0(temp.miss, ".*"))

  id1 <- data.frame(ID=miss$IID, F1=miss$F_MISS, stringsAsFactors = FALSE)
  id2 <- data.frame(ID=miss$IID, F2=miss$F_MISS, stringsAsFactors = FALSE)

  temp <- merge(temp, id2, by.x="ID2", by.y="ID", all.x = TRUE, sort=FALSE)
  temp <- merge(temp, id1, by.x="ID1", by.y = "ID", all.x = TRUE, sort = FALSE)

  #create second list resolving these pairs
  drop2 <- ifelse(temp$F1 > temp$F2, temp$ID1, temp$ID2)

  #create final list
  final <- unique(c(drop, drop2))
  print(paste("removing", length(final), "subjects resolving", nrow(rel_pairs), "relative pairs"))
  if(FID == TRUE){
    #get plink fam file in case there are family ids
    fam <- data.table::fread(paste0(plink.file, ".fam"), data.table = FALSE)
    fam <- fam[fam$V2 %in% final,]
    final <- fam[1:2]
  }

  return(final)
}
