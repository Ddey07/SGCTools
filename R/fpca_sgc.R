#' @title Functional principal component analysis of continuous, truncated and discrete functional data
#'
#' @description Covariance estimation and functional principal component analysis (FPCA) of continuous/truncated/ordinal/binary functional data based on the latent Semiparametric Gaussian copula process.
#'
#' @param X A numeric data matrix (n by m), n is the sample size and m is the number of time-points. Binary values must be coded as 0, 1. Ordinal values with l ordered categories must be coded as 0, 1, 2, ... , (l-1).
#' @param type A character denoting the type of variable in \code{X}, must be in the format "cont", "trunc", "ord", "bin".
#' @param argvals A vector of length m denoting the argument values of the functions. If NULL, we default it to an equidistant grid from 0 to 1.
#' @param df An integer denoting the degrees of freedom corresponding to the basis function.
#' @param T_out If supplied, the estimated covariance and FPCA is returned on this grid of argument values (time-points)
#' @param posd_tol If supplied, the posd_tol passed onto nearPD for nearest positive definite matrix
#' @param min_no_pairs The minimum number of pairs of observations used to calculate Kendall's Tau (otherwise report NA)
#' @param npc Prescribed value for the number of principal components. Defaults to 4.
#' @param scores A logical variable indicating if the user wants to return FPC scores or not.
#' @param impute A logical variable indicating whether missing values should be imputed before calculating FPC scores.
#' @param weights A logical variable indicating whether objective function needs to be weighted inversely proportional to sample size for each pair
#' @return \code{fpca.sgc.lat} returns
#' \itemize{
#'       \item{cov: }{the estimated m by m covariance matrix}
#'       \item{cov_out:}{If T_out is provided, the output covariance matrix calculated on the T_out grid}
#'       \item{efunctions: }{first npc number of eigenfunctions}
#'       \item{evalues: }{first npc number of eigenvalues}
#'       \item{latent: }{predictions of latent continuous trajectories, if \code{scores == TRUE}}
#'       \item{scores: }{first npc number of PC scores, if \code{scores == TRUE}}
#'       \item{par: }{Estimated parameters from the nls algorithm}
#'       \item{vcov:}{Variance-covariance matrix of the estimatd parameters from the nls algorithm}
#' }
#' @references
#' Dey D., Ghosal R., Merikangas K., Zipunnikov V. (2024) "Functional Principal Component Analysis for Continuous Non-Gaussian, Truncated, and Discrete Functional Data" Statistics in Medicine, 43, 5431-5445. \doi{10.1002/sim.10240}
#' @export
#' @import stats
#' @importFrom Matrix nearPD
#' @example man/examples/fpca_sgc_ex.R

fpca.sgc.lat = function(X, type,argvals=NULL, df = 5, T_out= NULL, npc = 4, scores = FALSE, impute=FALSE, min_no_pairs=30, weights = TRUE, posd_tol = 1e-08){

  # Check if X is a data frame or matrix
  if (!is.matrix(X)) {
    stop("Input 'X' must be a n * m matrix.")
  }

  # Get n (number of subjects) and m (number of time-points)
  n <- nrow(X)
  m <- ncol(X)

  # If argvals (length of time) is not supplied, generate it
  if (is.null(argvals)) {
    argvals = seq(0, 1, length = m)
  }

  # Check if df is a positive integer
  if (!is.numeric(df) || df <= 0 || df != round(df)) {
    stop("The 'df' parameter must be a positive integer.")
  }

  # Check if npc is a positive integer
  if (!is.numeric(npc) || npc <= 0 || npc != round(npc)) {
    stop("The 'npc' parameter must be a positive integer.")
  }

  # Check that 'type' is a single supported variable type
  if (length(type) != 1 || !type %in% c("cont", "bin", "ord", "trunc")) {
    stop("'type' must be one of \"cont\", \"bin\", \"ord\", \"trunc\".")
  }


  # calculate number of all pairs of time-points
  cmb <- combn(m,2)


  # create basis matrix
  knotsT = seq(min(argvals),max(argvals),l=df-2)
  norder = 4
  nbasisT = length(knotsT) + norder - 2
  dayrngT = c(min(argvals),max(argvals))
  bbasisT = create.bspline.basis(dayrngT,nbasisT,norder,knotsT)
  bs.argvals<-eval.basis(argvals,bbasisT)

  # calculate spline formula
  formula.spline <- function(j,l){paste(paste0("u",min(j,l),max(j,l)),"*",paste0("Tj",j),"*",paste0("Tl",l))}
  formula.spline <- Vectorize(formula.spline, vectorize.args = c("j","l"))
  spl.f <- paste(as.character(outer(1:df, 1:df, formula.spline)),collapse = ' + ')

  # Create data frame to solve the non-linear list squares problem
  if(type=="bin"){
    fjl.df <- function(j, l, data = X) {
      Tj = bs.argvals[j, ]
      Tl = bs.argvals[l, ]
      zratioj = mean(data[, j] == 0, na.rm = TRUE)
      zratiol = mean(data[, l] == 0, na.rm = TRUE)
      # count number of observations
      njl <- sum(!is.na(data[,j]*data[,l]))

      if(njl > min_no_pairs){

        tjl <- tryCatch(Kendall_mixed(cbind(data[,j],data[,l]))[1,2], error=function(e) NA)


        return(c(k=tjl, Tj, Tl, dj= zratioj, dl= zratiol, njl = njl))}
      else{
        return(c(k=NA, Tj, Tl, dj= zratioj, dl= zratiol, njl = njl))
      }
    }

    # fix bridging functions and create NLS formula
    bridgeF_bb_v <- Vectorize(bridgeF_bb,vectorize.args = c("r","zratio1","zratio2"))
    eunsc <- as.formula(paste0("k ~ ","bridgeF_bb_v(ginv(",spl.f, "), dj, dl)"))
    obj_df_colnames <- c("k",paste0("Tj",1:ncol(bs.argvals)), paste0("Tl",1:ncol(bs.argvals)),"dj", "dl","njl")

  } else if(type == "ord"){
    fjl.df <- function(j, l, data = X) {
      Tj = bs.argvals[j, ]
      Tl = bs.argvals[l, ]
      cats_j <- sort(unique(na.omit(data[, j])))
      hatdelta_j <- unlist(lapply(cats_j, function(x) {
        qnorm(1 - mean(data[, j] >= x, na.rm = TRUE))
      }))[-1]
      ncats_j = length(cats_j)
      cats_l <- sort(unique(na.omit(data[, l])))
      hatdelta_l <- unlist(lapply(cats_l, function(x) {
        qnorm(1 - mean(data[, l] >= x, na.rm = TRUE))
      }))[-1]
      ncats_l = length(cats_l)

      # count number of observations
      njl <- sum(!is.na(data[,j]*data[,l]))

      if(njl > min_no_pairs){

        tjl <- tryCatch(Kendall_mixed(cbind(data[,j],data[,l]))[1,2], error=function(e) NA)


        return(c(k=tjl, Tj, Tl, hatdelta_j, hatdelta_l, njl = njl))}
      else{
        return(c(k=NA, Tj, Tl, hatdelta_j, hatdelta_l, njl = njl))
      }
    }

    # fix bridging functions and create NLS formula
    ncutoff <- length(sort(unique(na.omit(X[,1]))))-1
    # defining bridging functions for our case
    l <- vector(mode="list",length =2*ncutoff + 1)
    names(l) = c("t",paste0("dj",1:ncutoff),paste0("dl",1:ncutoff))
    bridge_oo_fast_wrapper = make_function(args=l,body=quote({args_fa = unlist(as.list(environment())) # catch all the arguments as a list
    delta1 = args_fa[2:(ncutoff+1)] # assign values to baseline function arguments
    delta2 = args_fa[(ncutoff+2):(2*ncutoff + 1)]
    bridge_oo(args_fa[1],delta1,delta2)}))

    bridgeF_oo_v <- Vectorize(bridge_oo_fast_wrapper,vectorize.args = c("t",paste0("dj",1:ncutoff),paste0("dl",1:ncutoff)))
    eunsc <- as.formula(paste0("k ~ ","bridgeF_oo_v(ginv(",spl.f, "),",
                               paste(paste(paste0("dj",1:ncutoff),"=",paste0("dj",1:ncutoff),collapse=","),", ",
                                     paste(paste0("dl",1:ncutoff),"=",paste0("dl",1:ncutoff),collapse=",")),")"))
    obj_df_colnames <- c("k", paste0("Tj",1:ncol(bs.argvals)), paste0("Tl",1:ncol(bs.argvals)),paste0("dj",1:ncutoff), paste0("dl",1:ncutoff),"njl")

  } else if (type == "trunc"){
    fjl.df <- function(j, l, data = X) {
      Tj = bs.argvals[j, ]
      Tl = bs.argvals[l, ]
      zratioj = mean(data[, j] == 0, na.rm=TRUE)
      zratiol = mean(data[, l] == 0, na.rm=TRUE)
      # count number of observations
      njl <- sum(!is.na(data[,j]*data[,l]))

      if(njl > min_no_pairs){

        tjl <- tryCatch(Kendall_mixed(cbind(data[,j],data[,l]))[1,2], error=function(e) NA)


        return(c(k=tjl, Tj, Tl, dj= zratioj, dl= zratiol, njl = njl))}
      else{
        return(c(k=NA, Tj, Tl, dj= zratioj, dl= zratiol, njl = njl))
      }
    }

    # fix bridging functions and create NLS formula
    bridgeF_tt_v <- Vectorize(bridgeF_tt,vectorize.args = c("r","zratio1","zratio2"))
    eunsc <- as.formula(paste0("k ~ ","bridgeF_tt_v(ginv(",spl.f, "), dj, dl)"))
    obj_df_colnames <- c("k",paste0("Tj",1:ncol(bs.argvals)), paste0("Tl",1:ncol(bs.argvals)),"dj", "dl","njl")

  }

  # Get initial estimates assuming continuous data
  fjl.df.cont <- function(j, l, data = X) {
    Tj = bs.argvals[j, ]
    Tl = bs.argvals[l, ]
    # count number of observations
    njl <- sum(!is.na(data[,j]*data[,l]))

    if(njl > min_no_pairs){

      tjl <- tryCatch(Kendall_mixed(cbind(data[,j],data[,l]))[1,2], error=function(e) NA)


      return(c(k=tjl, Tj, Tl, njl = njl))}
    else{
      return(c(k=NA, Tj, Tl, njl = njl))
    }
  }

  # continuous bridging function
  bridgeF_cc_v <- Vectorize(bridgeF_cc,vectorize.args = c("r"))
  eunsc_cont <- as.formula(paste0("k ~ ","bridgeF_cc_v(ginv(",spl.f, "))"))
  obj_df_colnames_cont <- c("k",paste0("Tj",1:ncol(bs.argvals)), paste0("Tl",1:ncol(bs.argvals)), "njl")

  obj_df_cont <- data.frame(t(sapply(1:ncol(cmb), function(x){res = fjl.df.cont(cmb[1,x],cmb[2,x],data=X); if(is.na(res[1])){res = rep(NA,length(obj_df_colnames_cont))}; return(res)})))
  names(obj_df_cont) <- obj_df_colnames_cont

  # defining start point for NLS
  init <- runif(df*(df+1)/2,-5,5)
  formula.start <- function(j,l){paste(paste0("u",min(j,l),max(j,l)))}
  formula.start <- Vectorize(formula.start, vectorize.args = c("j","l"))
  names(init) <-  outer(1:df,1:df,formula.start)[lower.tri(outer(1:df,1:df,formula.start),diag=TRUE)]

  # run NLS optimization
  obj_df_cont = obj_df_cont[complete.cases(obj_df_cont),]
  if(weights==TRUE){
  ns0 <- nls(eunsc_cont,data=obj_df_cont,start=init, control=list(printEval=FALSE,warnOnly=TRUE, tol = 1e-03, minFactor = 1/32),weights= obj_df_cont$njl)
  } else{
    ns0 <- nls(eunsc_cont,data=obj_df_cont,start=init, control=list(printEval=FALSE,warnOnly=TRUE, tol = 1e-03, minFactor = 1/32))
  }
  ns1 <- ns0

  if (type != "cont"){
    # preparing data frame for NLS
    obj_df <- data.frame(t(sapply(1:ncol(cmb), function(x){res = fjl.df(cmb[1,x],cmb[2,x],data=X); if(is.na(res[1])){res = rep(NA,length(obj_df_colnames))}; return(res)})))
    names(obj_df) <- obj_df_colnames

    # defining start point for NLS
    init <- coef(summary(ns0))[,1]

    # run NLS optimization
    obj_df = obj_df[complete.cases(obj_df),]
    if(weights==TRUE){
    ns1 <- nlsLM(eunsc,data=obj_df,start=init, control=list(printEval=FALSE,warnOnly=TRUE, tol = 1e-03, minFactor = 1/32), trace = FALSE,  weights = obj_df$njl)
    } else {
      ns1 <- nlsLM(eunsc,data=obj_df,start=init, control=list(printEval=FALSE,warnOnly=TRUE, tol = 1e-03, minFactor = 1/32), trace = FALSE)
    }
  }

  uhat.nls <- matrix(ncol=df,nrow=df)
  uhat.nls[lower.tri(uhat.nls,diag=TRUE)] = coef(summary(ns1))[,1]
  uhat.nls[upper.tri(uhat.nls)] = t(uhat.nls)[upper.tri(uhat.nls)]

  # build correlation matrix on original scale
  Chat1 <- outer(argvals,argvals,Chat,u=uhat.nls, bs= bbasisT)
  Chat2 <- ginv(Chat1)
  diag(Chat2) <- 1
  Chat.grid0 <- nearPD(Chat2,corr=TRUE,maxit=10000,posd.tol = posd_tol)$mat
  Chat.grid <- Chat.grid0
  #Chat.grid.original <- Chat.grid2
  ee = eigen(Chat.grid)
  ee0 = ee
  res = list(cov = Chat.grid, efunctions = ee$vectors[,1:npc], evalues = ee$values[1:npc], vcov=vcov(ns1), par=ns1)


  #get smooth covariance matrix
  if(!is.null(T_out)){
    T1 = T_out
    knotsT = seq(min(T1),max(T1),l=df-2)
    norder = 4
    nbasisT = length(knotsT) + norder - 2
    dayrngT = c(min(T1),max(T1))
    bbasisT = create.bspline.basis(dayrngT,nbasisT,norder,knotsT)
    bs.argvals<-eval.basis(T1,bbasisT)

    Chat1 <- outer(T1,T1,Chat,u=uhat.nls,bs= bbasisT)
    Chat2 <- ginv(Chat1)
    diag(Chat2) <- 1
    Chat.grid <- nearPD(Chat2,corr=TRUE,maxit=1000,posd.tol = posd_tol)$mat
    #Chat.grid2 <- nearPD(Chat2,corr=TRUE,maxit=10000, posd.tol = 1e-03)$mat
    ee = eigen(Chat.grid)
    res = list(cov = as.matrix(Chat.grid0), cov_out = as.matrix(Chat.grid), efunctions = ee$vectors[,
                                                                                                    1:npc], evalues = ee$values[1:npc], vcov=vcov(ns1), par=ns1)
  }

  # check if data is dense or sparse (missing or not)
  dense = !any(is.na(X))


  # if we need to report PC scores
  if(scores==TRUE & impute == TRUE & dense == FALSE){
    L01 = getLatentPreds(X,type=rep(type,m), lat.cov.est = Chat.grid0, impute.missing = TRUE)
    res$latent = L01
    res$scores = L01 %*% ee0$vectors[,1:npc]
  } else if(scores == TRUE & dense == TRUE){
    L01 = getLatentPreds(X,type=rep(type,m), lat.cov.est = Chat.grid0, impute.missing = FALSE)
    res$latent = L01
    res$scores = L01 %*% ee0$vectors[,1:npc]
  }else if(scores==TRUE & impute == FALSE & dense == FALSE){
    L01 = getLatentPreds(X,type=rep(type,m), lat.cov.est = Chat.grid0, impute.missing = FALSE)
    res$latent = L01
    res$scores = t(sapply(1:n,function(i){sub= which(!is.na(X[i,])); cov_sub = Chat.grid0[sub,sub]; as.numeric(diag(ee0$values[1:npc]) %*% t(ee0$vectors[sub,1:npc]) %*% solve(cov_sub) %*% L01[i,sub])}))
  }
  return(res)
}
