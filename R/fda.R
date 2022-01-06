########################################
# Run FDA and get fitted smoothing curves
########################################

#usethis::use_package("fda")

#' Run FDA with a Fourier basis and positive constraints
#'
#' @param dat matrix of gene expressions with time x genes, size N x G
#' @param time vector of time
#' @param numbasis integer, number of Fourier basis
#' @param seed_num seed number
#' @return matrix of smoothing curves, size M x G
#' @export
#'
get_FDAcurves <- function(dat, time=NULL, numbasis=5, seed_num=12345){
  matrix_fda = as.matrix(dat)
  rownames(matrix_fda) = rownames(dat)
  if(is.null(time)){
    time = as.numeric(as.character(rownames(dat)))
  }
  #n_basis = 5

  if(max(time) - min(time) < 24){
    fda_time = seq(0, 23, by = 1)
  } else{
    fda_time = seq(min(time), max(time), by = 1)
  }

  set.seed(seed_num)
  # Fit a smooth line (FDA Fourier basis without positive constraints)
  basis = fda::create.fourier.basis(rangeval = range(fda_time),
                                    nbasis = numbasis,
                                    period=24)
  # define the harmonic acceleration differential operator
  harmLfd = fda::vec2Lfd(bwtvec=c(0, (2*pi/24)^2,0),
                         rangeval=c(0, max(fda_time)))

  # define the roughness penalty
  harmpar = fda::fdPar(fdobj=basis,
                       Lfdobj=harmLfd,
                       lambda = 1e-4)
  basic_basis = fda::smooth.basis(time, matrix_fda, harmpar)

  # Fit with positive constraints

  # redefine the roughness penalty
  harmpar = fda::fdPar(fdobj=basis,
                       Lfdobj=harmLfd,
                       lambda = 1e-10)

  # get coef of pos smoothing; since the default smooth.pos() function had issue with having multiple entries (multiple genes) at one time point
  pos.result.coef <- matrix(ncol = ncol(dat), nrow = numbasis)
  for (i in 1:ncol(dat)){
    data.to.pos = dat[,i]
    group.pos = fda::smooth.pos(time, data.to.pos, harmpar, dbglev=0) #suppress output=0
    pos.result.coef[,i] = group.pos$Wfdobj$coefs
  }
  # manipulate the pos.result$fd object obtained from smooth.pos() so it is in the same format as the fd obtained from smooth.basis
  group.pos$Wfdobj$coefs <- pos.result.coef
  group.pos$Wfdobj$fdnames$args <- basic_basis$fd$fdnames$time
  group.pos$Wfdobj$fdnames$reps <- basic_basis$fd$fdnames$reps
  group.pos$Wfdobj$fdnames$funs <- basic_basis$fd$fdnames$values

  group.pos.val = fda::eval.posfd(evalarg=fda_time, Wfdobj=group.pos$Wfdobj)
  colnames(group.pos.val) = colnames(dat)

  fda_curves = data.frame(fda_time, group.pos.val)

  return(fda_curves)
}
