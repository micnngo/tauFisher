
########################################
# Create difference matrix
########################################
#' Creates a matrix of gene expression  differences
#'
#' @param genes vector of genes of interest, length G
#' @param dat data matrix with time x genes, size N x G
#'
#' @importFrom magrittr %>%
#'
#' @return matrix of gene expression differences, size N x G^2
#' @export
#'
create_DiffMatrix <- function(genes, dat){

  combos <- expand.grid(genes, genes) %>%
    dplyr::mutate(same_genes = Var1==Var2) %>%
    dplyr::filter(same_genes == FALSE) %>%
    dplyr::mutate("Pair" = paste0(Var1, "_", Var2))
  diffs <- matrix(NA, nrow = nrow(combos), ncol = nrow(dat),
                  dimnames = list(combos$Pair, rownames(dat)))

  for(p in 1:nrow(combos)){
    for(t in 1:ncol(diffs)){
      diffs[p,t] <- dat[t, combos$Var1[p]] - dat[t, combos$Var2[p]]
    }
  }

  return(t(diffs))
}


########################################
# Filter out duplicate gene pairs
########################################

#' Filters difference matrix to keep only one pair of genes with differences > 0
#' Assumes effect of Gene1_Gene2 is not the same as Gene2_Gene1 on results
#'
#' @param diffs matrix of gene expression differences with time x gene pairs, size N x G^2
#' @return matrix of gene expression differences, size N x ((G^2 - G)/2)
#' @export
#'
filter_DiffMatrix <- function(diffs){

  # Filter out duplicate gene pairs
  same_genes <- apply(diffs, 2, sum)
  diffs <- diffs[, which(same_genes != 0)]

  # Filter out gene pairs with no differences
  # genes_min = apply(diffs, 2, min)
  # genes_max = apply(diffs, 2, max)
  # max_minus_min = genes_max - genes_min
  # diffs <- diffs[, which(max_minus_min != 0)]

  return(diffs)
}


########################################
# Scale matrix, row-wise [0,1]
########################################

#' Scales difference matrix row-wise (by time) to be between 0 and 1
#'
#' @param diffs matrix of gene expression differences with time x gene pairs, size N x P
#' @return matrix of gene expression differences, size N x P
#' @export
#'
scale_DiffMatrix <- function(diffs){

  diffs_scaled <- as.data.frame(t(apply(diffs, 1, scales::rescale)))

  return(diffs_scaled)
}
