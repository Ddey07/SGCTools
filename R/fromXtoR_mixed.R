#' @title Get an estimate of latent correlation matrix
#'
#' @description Estimation of latent correlation matrix from observed data of mixed types (continuous/truncated/ordinal/binary) based on the Semiparametric Gaussian copula model.
#'
#' @param X A numeric data matrix (n by p), n is the sample size and p is the number of variables. Binary variables must be coded as 0, 1. Ordinal variables with l ordered categories must be coded as 0, 1, 2, ... , (l-1).
#' @param type A vector of length p containing type of variables in \code{X}, must be in the format "cont", "trunc", "ord", "bin". If not supplied, our function tries to guess the data-type from the observations
#' @param use.nearPD A logical value indicating whether to use \link[Matrix]{nearPD} or not when the resulting latent correlation matrix is not positive definite.
#' @param deriv A logical value indicating whether to calculate point-wise derivatives of the bridging functions useful for calculating the asymptotic confidence interval of latent correlations
#' @return \code{fromXtoRmixed} returns
#' \describe{
#'       \item{hatR}{the estimated latent correlation matrix}
#'       \item{hatRprime}{if (deriv = TRUE), returns the pointwise derivative of the bridging function with respect to latent correlation}
#' }
#' @references
#' Dey D., Zipunnikov V. (2022) "Semiparametric Gaussian Copula Regression modeling for Mixed Data Types (SGCRM)" <https://arxiv.org/abs/2205.06868>
#'
#' Fan J., Liu H., Ning Y. and Zou H. (2017) "High dimensional semiparametric latent graphicalmodel for mixed data" <doi:10.1111/rssb.12168>.
#'
#' Yoon G., Carroll R.J. and Gaynanova I. (2020) "Sparse semiparametric canonical correlation analysis for data of mixed types" <doi:10.1093/biomet/asaa007>.
#'
#' Yoon G., Mueller C.L., Gaynanova I. (2020) "Fast computation of latent correlations" <https://arxiv.org/abs/2006.13875>.
#'
#' @export
#' @import stats
#' @importFrom Matrix nearPD
#' @example man/examples/fromXtoR_ex.R
fromXtoRMixed=function(X,type=NULL,use.nearPD= TRUE,deriv=FALSE){
  # Check if X is a matrix or data frame
  if (!is.data.frame(X) && !is.matrix(X)) {
    stop("X must be a data frame or a matrix.")
  }

  n = nrow(X)
  p = ncol(X)

  #Classifying column types
  if(is.null(type)){
    bin_id=which(apply(X,2,function(x) { all(na.omit(x) %in% 0:1) }))
    ord_id=which(apply(X,2,function(x) { sum(na.omit(x) %%1)==0 && length(unique(na.omit(x[x!=0]))) >=2 && length(unique(na.omit(x[x!=0]))) <15}))
    trunc_id=which(apply(X,2,function(x) { x=na.omit(x); mean(x==0) > 0.01 && sum(x %% 1) > 0 }))
    type <- rep("cont",p)
    type[bin_id] <- "bin"
    type[ord_id] <- "ord"
    type[trunc_id] <- "trunc"
  }

  # pre-process data to estimate cutoffs and transformation functions
  PR=preprocess_data(X)$features

  hatR=matrix(1,p,p)
  hatRprime <- hatR

    if(deriv==TRUE){
      for(i in 1:(p-1)){
        for(j in (i+1):p){
          pr0 <- PR[c(i,j)]
          h = fromXtoR_bivariate(X[,i],X[,j],PR=pr0,type1=type[i],type2=type[j],deriv=deriv)
          hatR[i,j]=h[1]
          hatRprime[i,j]=h[2]
          hatR[j,i]=hatR[i,j]
          hatRprime[j,i]=hatRprime[i,j]
        }
      }
    } else {
      for(i in 1:(p-1)){
        for(j in (i+1):p){
          pr0 <- PR[c(i,j)]
          h = fromXtoR_bivariate(X[,i],X[,j],PR=pr0,type1=type[i],type2=type[j],deriv=deriv)
          hatR[i,j]=h[1]
          hatR[j,i]=hatR[i,j]
        }
      }
    }

  if(use.nearPD==TRUE){
    hatR=nearPD(hatR, corr=TRUE, maxit=1000, posd.tol = 1e-3)$mat
  }

  # convert hatR to matrix object
  hatR= as.matrix(hatR)

  if(deriv==TRUE){
    return(list(hatR=hatR,hatRprime=hatRprime))
  } else {
    return(list(hatR=hatR))
  }
}
