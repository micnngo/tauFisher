

########################################
# Run MetaCycle
########################################


#' Runs the meta2d function from MetaCycle
#'
#' @param df expression data matrix with genes as rows and time/samples as columns
#' @param timepoints vector of time points, equal size to number of columns of df
#' @param out_file name of file in ".txt" format; saves expression matrix into text file for input into MetaCycle
#' @param out_dir output directory to save MetaCycle results
#' @param method character vector defining which method(s) to detect rhythmic signals
#' @param ... additional parameters for meta2d from MetaCycle
#' @export
#'
run_metacycle <- function(df,
                          timepoints,
                          out_file,
                          out_dir,
                          method=c("JTK", "LS"),
                          ...){

  # Need to write data files into text files for MetaCycle
  dfn = cbind("Genes"=rownames(df), df)
  utils::write.table(dfn, file=out_file,sep="\t", quote=FALSE, row.names=FALSE)

  MetaCycle::meta2d(infile=out_file,
                    filestyle="txt",
                    outdir=out_dir,
                    timepoints=timepoints,
                    minper = 20, maxper = 28,
                    cycMethod = method,
                    outIntegration = "both",
                    ...)
  cat("Finished finding rhythmic genes!")
}



########################################
# Get the periodic genes
########################################

#' Extracts the significant periodic genes from the meta2d output
#' @param input_file input file to read in, needs to be ".txt" (output of meta2d)
#' @param test_genes vector of genes in the testing data set
#' @param per_method vector of periodic method to use (LS or JTK)
#' @param thres threshold for p-value
#' @importFrom magrittr %>%
#' @export
#'
find_periodic_genes <- function(input_file,
                                test_genes,
                                per_method=c("LS", "JTK"),
                                thres=1.01){

  cc_genes = toupper(c("Bmal1", "Arntl", "Dbp", "Nr1d1", "Nr1d2","Per1", "Per2", "Per3", "Cry1",  "Cry2"))

  test_genes = toupper(test_genes)

  if(any(per_method %in% "LS")){

    per_LS <- utils::read.delim(input_file) %>%
      dplyr::filter(round(LS_period) == 24, LS_pvalue < thres) %>%
      dplyr::arrange(LS_pvalue) %>%
      dplyr::mutate(CycID = toupper(CycID))

    genes_LS = per_LS[per_LS$CycID %in% test_genes, ] %>%
      dplyr::top_n(-10, LS_pvalue) %>%
      dplyr::select(CycID) %>% c(.)
    genes_LS = unique(c(genes_LS$CycID, cc_genes[cc_genes %in% per_LS$CycID]) )
    #genes_LS = unique(c(genes_LS$CycID, cc_genes) )


  } else{ genes_LS = NULL}

  if(any(per_method %in% "JTK")){

    per_JTK <- utils::read.delim(input_file) %>%
      dplyr::filter(round(JTK_period) == 24, JTK_pvalue < thres) %>%
      dplyr::arrange(JTK_pvalue) %>%
      dplyr::mutate(CycID = toupper(CycID))

    genes_JTK = per_JTK[per_JTK$CycID %in% test_genes, ] %>%
      dplyr::top_n(-10, JTK_pvalue) %>%
      dplyr::select(CycID) %>% c(.)
    genes_JTK = unique(c(genes_JTK$CycID, cc_genes[cc_genes %in% per_JTK$CycID]))
    #genes_JTK = unique(c(genes_JTK$CycID, cc_genes))



  } else{ genes_JTK = NULL}


  if(any(!(per_method %in% c("JTK", "LS")))){
    print("Not a valid method included. Must be at least one of JTK or LS.")
  }

  return(list(LS = genes_LS, JTK = genes_JTK))
}

