#' @title Get predictors of latent variables from observed data
#'
#' @description BLUPs of latent variables from observed data and estimated latent correlation matrix
#' @param X A numeric data matrix (n by p), n is the sample size and p is the number of variables. Binary variables must be coded as 0,1. Ordinal variables with l ordered categories must be coded as 0, 1, 2, ... , (l-1).
#' @param type A vector of length p containing type of variables in \code{X}, must be in the format "cont", "trunc", "ord", "bin". If not supplied, our function tries to guess the data-type from the observations.
#' @param lat.cov.est The estimate of the latent correlation matrix from \code{fromXtoRMixed}
#' @param impute.missing A logical variable indicating whether the missing data will be imputed with Missing at Random assumptions or not.
#' @param ridge A non-negative scalar added to the diagonal of the observed-pattern submatrices before each per-subject BLUP solve. This regularizes near-singular solves that can otherwise inflate individual latent predictions to extreme values when the estimated latent correlation is poorly conditioned (e.g. for sparse or binary margins). Defaults to 0.05; set to 0 to reproduce the unregularized behaviour.
#' @param pmv_maxpts,pmv_abseps Monte-Carlo controls for the truncated-normal conditional mean of the observed binary/ordinal cells (\code{tmvtnorm::mtmvnorm} via \code{mvtnorm::GenzBretz}). Defaults (25000, 1e-3) reproduce the original precision. For large-scale latent prediction the conditional mean is averaged over many cells, so a much cheaper setting (e.g. \code{pmv_maxpts = 500}, \code{pmv_abseps = 1e-2}) gives a large speed-up at negligible accuracy cost.
#' @return \code{getLatentPreds} returns a n by p matrix containing continuous predictions of observed mixed latent variables.
#' @references
#' Dey ., Zipunnikov V. (2022) "Semiparametric Gaussian Copula Regression modeling for Mixed Data Types (SGCRM)" <https://arxiv.org/abs/2205.06868>
#' @export
#' @import stats
#' @importFrom Matrix nearPD
#' @example man/examples/fromXtoR_ex.R
getLatentPreds = function(X, type = NULL, lat.cov.est, impute.missing = TRUE, ridge = 0.05,
                          pmv_maxpts = 25000L, pmv_abseps = 1e-3){
  PRx = preprocess_data(X)
  PR=PRx$features
  L0 = X
  for(i in 1:nrow(X)){
    L0[i,]=recover_row(i,X=X,lat.cov = lat.cov.est, feature=PR, type = type, impute = impute.missing,
                       ridge = ridge, pmv_maxpts = pmv_maxpts, pmv_abseps = pmv_abseps)
  }
  return(L0)
}

