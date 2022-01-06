
########################################
# Train tauFisher
########################################
#' Train tauFisher
#'
#' @param periodic_genes vector of periodic genes to subset on
#' @param train_df gene expression matrix for training data (gene x sample/time)
#' @param train_time vector of time points (e.g. c(0, 4, 8, ..., 48, 52))
#' @param nrep number of replicates
#' @param numbasis number of basis functions for FDA
#'
#' @importFrom magrittr %>%
#'
#' @return list of (1) periodic genes in the training set,
#'                 (2) PCA embedding,
#'                 (3) regression model
#' @export
#'
train_tauFisher <- function(periodic_genes,
                            train_df,
                            train_time,
                            nrep=1,
                            numbasis=5){

  # Subset data on periodic genes
  rownames(train_df) = toupper(rownames(train_df))
  chosen_genes <- periodic_genes[periodic_genes %in% rownames(train_df)]

  df <- t(data.frame(train_df[chosen_genes, ]))
  rownames(df) <- train_time

  # Order the genes
  chosen_genes <- chosen_genes[order(chosen_genes)]
  df = df[, order(colnames(df))]

  # log2 transform
  df_log = log2(df+1)

  # Convert time to 24 hours
  #train_time24 = convert_to_24hr(train_time)

  # Run FDA
  fda_expr = get_FDAcurves(df_log, time=train_time, numbasis=numbasis)  %>%
    #dplyr::filter(!(fda_time %% 1)) %>%
    dplyr::mutate(time_24 = fda_time - 24*floor(fda_time/24)) %>%
    dplyr::mutate(time_lab = paste0(time_24, "_",
                                    rep((max(nrep)+1):(max(nrep)+max(nrep)),
                                        each = 24, length = nrow(.))))
  new_fda_rownames = fda_expr$time_lab

  # Calculate differences
  fda_expr2 = fda_expr[, -c(1, ncol(fda_expr)-1, ncol(fda_expr))]
  #print(fda_expr2[1:5,1:5])
  fda_diff <- create_DiffMatrix(chosen_genes, fda_expr2)
  #fda_diff_filtered = filter_DiffMatrix(fda_diff)
  fda_diff_scaled = scale_DiffMatrix(fda_diff)

  fda_mat = as.matrix(fda_diff_scaled)

  # Train data
  train = fda_mat
  rownames(train) = new_fda_rownames
  train_time = fda_expr$time_24

  # PCA
  X_PCA <- train
  X_PCA <- as.data.frame(X_PCA)
  X_PCA$CT <- train_time

  pc <- stats::prcomp(X_PCA[, -ncol(X_PCA)], scale = FALSE)
  #pc_pred <- predict(pc, newdata = test)

  # Multinomial regression
  ndims = 2
  pc_data <- data.frame(pc$x[,1:ndims])
  #pc_data$CT24 <- as.numeric(vapply(stringr::str_split(rownames(pc_data), pattern = '_'),
  #                                  '[', 1, FUN.VALUE = character(1) ))
  pc_data$CT24 <- train_time
  pc_data$CT24_rl <- stats::relevel(factor(pc_data$CT24),
                                    ref = as.character(min(train_time)))

  mod <- nnet::multinom(CT24_rl ~ PC1 + PC2, data = pc_data, trace=T)
  #pred_vals = stats::predict(m1, newdata = pc_pred[,1:ndims])
  #pred_vals = as.numeric(as.character(pred_vals))

  #out_df = data.frame("truth" = test_time, "predicted" = pred_vals)

  # just in case a gene gets filtered out (not likely to occur)
  # usually, filtered_chosen_genes = chosen_genes
  #gene_pairs = colnames(fda_diff_filtered)
  #filtered_chosen_genes = sort(unique(unlist(stringr::str_split(gene_pairs, "_"))))

  return(list(Genes=chosen_genes, PCA_embedding=pc, Model=mod))
}

########################################
# Test tauFisher
########################################

#' Test tauFisher
#'
#' @param chosen_genes vector of chosen periodic genes to subset on
#' @param test_df gene expression matrix for training data (gene x sample/time)
#' @param test_time vector of times (e.g. c(0, 2, 4, ..., 26))
#' @param pca_embedding matrix of PCs and loadings from `train_tauFisher`
#' @param trained_model trained model from `train_tauFisher`
#'
#' @importFrom magrittr %>%
#'
#' @return data frame with true time and predicted time
#' @export
#'
test_tauFisher <- function(chosen_genes,
                           test_df,
                           test_time,
                           pca_embedding,
                           trained_model){

  # Subset data on chosen genes (periodic genes from train set)
  rownames(test_df) = toupper(rownames(test_df))
  df <- t(data.frame(test_df[chosen_genes, ]))

  # Order the genes
  chosen_genes <- chosen_genes[order(chosen_genes)]
  df = df[, order(colnames(df))]

  # log2 transform
  df_log = log2(df+1)

  # Test data
  # Create differences matrix
  df_diff <- create_DiffMatrix(chosen_genes, df_log)
  #df_diff_filtered = filter_DiffMatrix(df_diff)
  df_diff_scaled = scale_DiffMatrix(df_diff)

  test = df_diff_scaled

  # PCA
  pc_pred <- stats::predict(pca_embedding, newdata = test)

  # Multinomial regression
  ndims = 2
  pred_vals = stats::predict(trained_model, newdata = pc_pred[,1:ndims])
  pred_vals = as.numeric(as.character(pred_vals))

  out_df = data.frame("truth" = test_time, "predicted" = pred_vals)

  return(out_df)
}
