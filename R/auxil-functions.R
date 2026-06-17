# import function from namespace
fromKtoR_mixed = getFromNamespace("fromKtoR_mixed", "mixedCCA")

# Make function using argument and body
make_function <- function(args, body, env = parent.frame()) {
  args <- as.pairlist(args)

  eval(call("function", args, body), env)
}

# transformation function: inverse of the Fisher-style map g(x) = log((1+x)/(1-x));
# maps the spline linear predictor onto the correlation scale (-1, 1).
ginv = function(x){
  (exp(x)-1)/(exp(x)+1)
}


# Calculate estimated covariance matix for spline
Chat <- function(u,s,t,bs){
  fb.sp(u, as.numeric(eval.basis(s,bs)), as.numeric(eval.basis(t,bs)))
}
Chat <- Vectorize(Chat,vectorize.args = c("s","t"))


# calculate concordances
conc_ties <- function(X){
  if (!is.matrix(X) && !is.data.frame(X)) stop ("X must be a matrix or data.frame.")
  if (!nrow(X) > 1 || !ncol(X) > 1) stop ("X must have at least 2 rows and 2 columns.")

  n <- as.integer(nrow(X)); d <- ncol(X)
  C <- D <- array(NA, dim = c(n,d,d))
  N <- diag(d)

  for (i in 1:(d-1)){
    keep1 <- which(!is.na(X[,i]))
    ord <- order(X[keep1, i])
    for (j in (i+1):d){
      keep2_l <- !is.na(X[keep1,j])
      keep2_o <- which(keep2_l[ord])
      N[i, j] <- N[j, i] <- length(keep2_o)
      C[keep1[keep2_o], i, j] <- C[keep1[keep2_o], j, i] <- countSwaps(-X[keep1[ord],j][keep2_o], X[keep1[ord],i][keep2_o])
      D[keep1[keep2_o],i,j] <- D[keep1[keep2_o], j, i] <- countSwaps(X[keep1[ord],j][keep2_o], X[keep1[ord],i][keep2_o])
    }
  }
  return(list(C=C, D=D, N=N))
}

# NOTE: Kendall_mixed() is defined and exported in R/Kendall_mixed.R (the canonical
# version, which sets diag(tau) <- 1 and validates input). It is intentionally not
# redefined here; an earlier duplicate in this file shadowed the exported one.

# loading ord and truncated marginal bridging functions from mixedCCA
bridgeF_tt = getFromNamespace("bridgeF_tt", "mixedCCA")
bridgeF_cc = function(r){2/pi * asin(r)}

bridgeF_bb <- function (r, zratio1, zratio2)
{
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)
  if(r[1]==1){r=r-1e-6}
  if(r[1]==(-1)){r=r+1e-6}
  res <- as.numeric(2 * (fMultivar::pnorm2d(de1, de2, rho = r) -
                           zratio1 * zratio2))
  return(res)
}

# continuous - ord bridging function
bridge_co = function(t,delta,...){
  M <- 12
  delta <- c(delta, M)
  S <- cbind(c(1,0,t/sqrt(2)),c(0,1,-t/sqrt(2)), c(t/sqrt(2),-t/sqrt(2),1))
  l <- length(delta)
  term1 <- sum(sapply(1:(l-1), function(x){4*pmvnorm(upper=c(delta[x],delta[x+1],0),sigma=S,algorithm = Miwa()) - 2*pnorm(delta[x])*pnorm(delta[x+1])}))
  return(term1)
}


# ord - ord bridging function
bridge_oo = function(t,delta1,delta2){
  M <- 12
  l1 <- length(delta1)
  l2 <- length(delta2)
  if(t[1]==1){t=t-1e-6}
  if(t[1]==(-1)){t=t+1e-6}
  setting1 <- expand.grid(c(1:l1),c(1:l2))
  delta_t1 <- c(-M,delta1,M)
  delta_t2 <- c(-M,delta2,M)
  #delta_t3 <- c(-M,delta2)

  term1= sum(sapply(1:nrow(setting1), function(i){x <- as.numeric(setting1[i,]) ;
  pnorm2d(delta_t1[x[1]+1],delta_t2[x[2]+1], rho=t)*
    (pnorm2d(delta_t1[x[1]+2],delta_t2[x[2]+2],rho=t) -  pnorm2d(delta_t1[x[1]+2],delta_t2[x[2]],rho=t))}))

  term2 = sum(sapply(1:l1, function(i){x <- as.numeric(setting1[i,]) ;
  pnorm(delta_t1[x[1]+1])*(pnorm2d(delta_t1[x[1]+2],delta_t2[l2+1],rho=t))
  }))
  return(c(2*(term1-term2)))
}

# bin-ord bridging function
bridge_bo = function(t,delta1,delta2){
  M <- 12
  if(t[1]==1){t=t-1e-6}
  if(t[1]==(-1)){t=t+1e-6}
  l1 <- length(delta1)
  delta_t1 <- c(delta2,M)
  delta_t2 <- c(-M,delta2)
  Splus <- cbind(c(1,t),c(t,1))
  Sminus <- cbind(c(1,-t),c(-t,1))

  term1 <- sum(sapply(1:(length(delta_t1)-1), function(x){pmvnorm(lower=c(delta_t1[x],-M),upper=c(delta_t1[x+1],-delta1),sigma=Sminus,algorithm = Miwa())*
      pmvnorm(upper=c(delta_t1[x],delta1),sigma=Splus,algorithm = Miwa())}))
  #term2 <- sum(sapply(1:(length(delta_t2)-1), function(x){pmvnorm(lower=c(delta_t2[x],-M),upper=c(delta_t2[x+1],-delta1),sigma=Sminus)*pmvnorm(upper=c(-delta_t2[x+1],delta1),sigma=Sminus)}))
  term3 <- sum(sapply(1:(length(delta_t1)-1), function(x){pmvnorm(lower=c(delta_t1[x],-M),upper=c(delta_t1[x+1],delta1),sigma=Splus,algorithm = Miwa())*pmvnorm(upper=c(delta_t1[x],-delta1),sigma=Sminus,algorithm = Miwa())}))
  return(c(2*(term1-term3)))
}


# ord - truncated bridging function
bridge_ot=function(t,delta1,delta2){
  M <- 12
  l1 <- length(delta1)
  delta_t1 <- c(delta1,M)
  delta_t2 <- c(-M,delta1)
  Splus <- cbind(c(1,t),c(t,1))
  Sminus <- cbind(c(1,-t),c(-t,1))
  S4plus <- cbind(c(1,0,0,-t/sqrt(2)),c(0,1,-t,t/sqrt(2)),c(0,-t,1,-1/sqrt(2)),c(-t/sqrt(2),t/sqrt(2),-1/sqrt(2),1))
  S4minus <- cbind(c(1,0,0,-t/sqrt(2)),c(0,1,t,-t/sqrt(2)),c(0,t,1,-1/sqrt(2)),c(-t/sqrt(2),-t/sqrt(2),-1/sqrt(2),1))

  term1 <- sum(sapply(1:(length(delta_t1)-1), function(x){pmvnorm(lower=c(delta_t1[x],-M),upper=c(delta_t1[x+1],-delta2),sigma=Sminus,algorithm = Miwa())*pmvnorm(upper=c(delta_t1[x],delta2),sigma=Splus,algorithm = Miwa()) +
      pmvnorm(lower=c(delta_t1[x],-M,-M,-M),upper=c(delta_t1[x+1],delta_t1[x],-delta2,0),sigma=S4plus,algorithm = Miwa())}))

  term2 <- sum(sapply(1:(length(delta_t2)-1), function(x){pmvnorm(lower=c(delta_t2[x],-M),upper=c(delta_t2[x+1],-delta2),sigma=Sminus,algorithm = Miwa())*pmvnorm(upper=c(-delta_t2[x+1],delta2),sigma=Sminus,algorithm = Miwa()) +
      pmvnorm(lower=c(delta_t2[x],-M,-M,-M),upper=c(delta_t2[x+1],-delta_t2[x+1],-delta2,0),sigma=S4minus,algorithm = Miwa())}))

  return(2*(term1-term2))
}


# Observed data to latent correlation estimation for two variables
fromXtoR_bivariate=function(X1,X2,PR,type1="ord",deriv=TRUE,type2="ord",tol=1e-6){

  ### Lexicographically ordering variable type for relatively easy program implementation
  type=c(type1,type2)
  X=cbind(X1,X2)[,order(type)]
  PR= PR[order(type)]
  type1=sort(type)[1]
  type2=sort(type)[2]
  K <- Kendall_mixed(X)[1,2]


  if(type1=="ord" & type2=="ord"){
      hatdelta1 <- PR[[1]]$cutoff[-c(1,length(PR[[1]]$cutoff))]
      hatdelta2 <- PR[[2]]$cutoff[-c(1,length(PR[[2]]$cutoff))]
      fitf_oo <- function(t){(bridge_oo(t,hatdelta1,hatdelta2)-K)^2}
      R=optimize(fitf_oo,c(-1,1))$minimum
      if(deriv==TRUE){
        Rprime <- 1/grad(function(t){bridge_oo(t,hatdelta1,hatdelta2)},R,method="simple")
      }
  }
  else if((type1=="ord" & type2=="trunc") | (type2=="ord" & type1=="trunc")){

      hatdelta1 <- PR[[1]]$cutoff[-c(1,length(PR[[1]]$cutoff))]
      hatdelta2 <- PR[[2]]$cutoff[-1]
      fitf_ot <- function(t){(bridge_ot(t,hatdelta1,hatdelta2)-K)^2}
      R=optimize(fitf_ot,c(-1,1))$minimum
      if(deriv==TRUE){
        Rprime <- 1/grad(function(t){bridge_ot(t,hatdelta1,hatdelta2)},R,method="simple")
      }

  }
  else if((type1=="cont" & type2=="ord") | (type2=="cont" & type1=="ord")){
    #cats <- sort(unique(na.omit(X[,2])))
      hatdelta1 <- PR[[2]]$cutoff[-c(1,length(PR[[2]]$cutoff))]
      fitf_co <- function(t){(bridge_co(t,hatdelta1)-K)^2}
      R=optimize(fitf_co,c(-1,1))$minimum
      if(deriv==TRUE){
        Rprime <- 1/grad(function(t){bridge_co(t,hatdelta1)},R,method="simple")
      }
  }
  else if((type1=="bin" & type2=="ord") |(type2=="bin" & type1=="ord")){
      hatdelta1 <- PR[[1]]$cutoff[-c(1,length(PR[[1]]$cutoff))]
      hatdelta2 <- PR[[2]]$cutoff[-c(1,length(PR[[2]]$cutoff))]
      fitf_bo <- function(t){(bridge_bo(t,hatdelta1,hatdelta2)-K)^2}
      R=optimize(fitf_bo,c(-1,1))$minimum
      if(deriv==TRUE){
        Rprime <- 1/grad(function(t){bridge_bo(t,hatdelta1,hatdelta2)},R,method="simple")
      }
  }
  else if(type1 %in% c("bin","cont","trunc") & type2 %in% c("bin","cont","trunc")){
    zratio1=NULL
    zratio2=NULL
    if(mean(X[,1]==0, na.rm=TRUE) > 0.02){zratio1=mean(X[,1]==0, na.rm=TRUE)}
    if(mean(X[,2]==0, na.rm=TRUE) > 0.02){zratio2=mean(X[,2]==0, na.rm=TRUE)}
    if(type1=="cont"){type1="continuous"}
    if(type2=="cont"){type2="continuous"}
    if(type1=="bin"){type1="binary"}
    if(type2=="bin"){type2="binary"}
    R=fromKtoR_mixed(K,type1=type1,type2=type2, zratio1=zratio1,zratio2=zratio2)
    fitf <- function(x){fromKtoR_mixed(x,type1=type1,type2=type2, zratio1=zratio1,zratio2=zratio2)}
    if(deriv==TRUE){
      Rprime <- grad(fitf,K,method="simple")
    }
  }
  if(deriv==TRUE){
    return(as.numeric(c(R,Rprime)))
  } else {
    return(R)
  }
}


npn.f.est <- function(x){
  x <- na.omit(x)
  En <- ecdf(x)
  n <- length(x)
  Dn <- 0
  # Dn <- 1/(4*n^(1/4)*sqrt(pi * log(n)))
  Fn <- function(x){
    if(En(x) > Dn & En(x) < 1-Dn){
      En(x)
    }
    else if(En(x) <= Dn){
      Dn
    }
    else{
      1-Dn
    }
  }

  mu <- mean(x)
  sigma <- sqrt(var(x))


  f.hat <- function(x){
    qnorm((n/(n+1)*Fn(x)))
  }
  return(f.hat)
}

# Preprocess observed data - calculate cutoffs, estimate transformation functions
preprocess_data=function(X,type=NULL,...){
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
  type1 <- type
  features <- list()
  for(i in 1:ncol(X)){
    if(type[i]=="cont"){
      features[[i]]=list(fn=npn.f.est(X[,i]))
    } else if(type[i]=="trunc"){
      features[[i]]=list(fn=npn.f.est(X[,i]),cutoff=c(-Inf,qnorm(mean(X[,i]==0,na.rm=TRUE))))
    } else if(type[i]=="bin"){
      features[[i]]=list(cutoff=c(-Inf,qnorm(1-mean(X[,i],na.rm=TRUE)),Inf))
    } else if(type[i]=="ord"){
      cats <- sort(unique(na.omit(X[,i])))
      hatdelta1 <- unlist(lapply(cats,function(x){qnorm(1-mean(X[,i]>=x, na.rm=TRUE))}))[-1]
      features[[i]]=list(cutoff=c(-Inf,hatdelta1,Inf))
    }
  }
  return(list(type=type1,features=features))
}

# recover latent data from observed data
recover_row=function(i,X=Z1,lat.cov,feature=NULL,type=NULL,impute=TRUE,ridge=0.05,...){
  n = nrow(X)
  p = ncol(X)

  if(is.null(feature)){
    feature=preprocess_data(X)
  }

  R2=lat.cov
  # ps stability fix: ridge the latent correlation before it feeds the
  # truncated-normal conditional moments (mtmvnorm) and the BLUP/imputation
  # solves. A poorly-conditioned, low-rank correlation (typical of the least-
  # informative margins, e.g. a binary component) otherwise yields extreme
  # conditional means for individual subjects (|z| up to ~1e2-1e69), which then
  # dominate cov(z) and collapse the downstream pfpca rank. cov2cor keeps the
  # latent variances at 1 so the thresholds/cutoffs stay valid. ridge = 0
  # reproduces the exact original behaviour.
  if(ridge > 0) R2 = stats::cov2cor(as.matrix(R2) + ridge * diag(nrow(as.matrix(R2))))

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


  M=which(is.na(X[i,]))
  O=which(!is.na(X[i,]))

  O_B= O[which(O %in% which(type %in% c("bin" , "ord")))]
  O_C= O[which(O %in% which(type == "cont"))]
  O_T=O[which(O %in% which(type == "trunc"))]


  Z=rep(NA,p)

  if(length(O_T)!=0){

    O_B=c(O_B,O_T[which(X[i,O_T]==0)])
    O_C=c(O_C,O_T[which(X[i,O_T]!=0)])

    Z[O_C]=as.numeric(lapply(O_C,function(x){f=feature[[x]]$fn; f(X[i,x])}))
    Z1 <- Z

    if(length(O_B)!=0){
      if(length(O_C)!=0){
        cond.mean = as.numeric((R2[O_B,O_C]) %*% solve(R2[O_C,O_C]) %*% (Z[O_C]))
        cond.sigma= R2[O_B,O_B] - R2[O_B,O_C] %*% solve(R2[O_C,O_C]) %*% R2[O_C,O_B]
      } else {
        cond.mean= rep(0,length(O_B))
        cond.sigma= R2[O_B,O_B]
      }

      lb=rep(-Inf,length(O_B))
      ub=rep(Inf,length(O_B))

      for(j in 1:length(O_B)){
        lb[j]=feature[[O_B[j]]]$cutoff[X[i,O_B[j]]+1]
        ub[j]=feature[[O_B[j]]]$cutoff[X[i,O_B[j]]+2]
      }


      Z1[O_B]=mtmvnorm(sigma = cond.sigma, mean=cond.mean,
                       lower = lb, upper=ub,
                       doComputeVariance=FALSE,
                       pmvnorm.algorithm=GenzBretz()) $tmean
    }
    #Imputing data
    if(impute==TRUE){
      if(length(M)>0){
        Z1[M] = as.numeric(R2[M,c(O_B,O_C)] %*% solve(R2[c(O_B,O_C),c(O_B,O_C)]) %*% (Z1[c(O_B,O_C)]))
      }
    }

  } else {
    Z[O_C]=as.numeric(lapply(O_C,function(x){f=feature[[x]]$fn; f(X[i,x])}))
    Z1 <- Z
    if(length(O_B)!=0){
      if(length(O_C)!=0){
        cond.mean = as.numeric((R2[O_B,O_C]) %*% solve(R2[O_C,O_C]) %*% (Z[O_C]))
        cond.sigma= R2[O_B,O_B] - R2[O_B,O_C] %*% solve(R2[O_C,O_C]) %*% R2[O_C,O_B]
      } else {
        cond.mean= rep(0,length(O_B))
        cond.sigma= R2[O_B,O_B]
      }

      lb=rep(-Inf,length(O_B))
      ub=rep(Inf,length(O_B))

      for(j in 1:length(O_B)){
        lb[j]=feature[[O_B[j]]]$cutoff[X[i,O_B[j]]+1]
        ub[j]=feature[[O_B[j]]]$cutoff[X[i,O_B[j]]+2]
      }


      Z1[O_B]=mtmvnorm(sigma = cond.sigma, mean=cond.mean,
                       lower = lb, upper=ub,
                       doComputeVariance=FALSE,
                       pmvnorm.algorithm=GenzBretz()) $tmean
    }
    #Imputing data
    if(impute==TRUE){
      if(length(M)>0){
        Z1[M] = as.numeric(R2[M,c(O_B,O_C)] %*% solve(R2[c(O_B,O_C),c(O_B,O_C)]) %*% (Z1[c(O_B,O_C)]))
      }
    }
  }
  return(Z1)
}

### covariance smooth function ###
# tensor-product spline evaluation used by Chat()
fb.sp <- function(u,Tj,Tl){
  t(Tj) %*% u %*% Tl
}


# asymptotic variance functions

lo.elim <- function(p){
  E_l <- matrix(0,nrow=p*(p-1)/2,ncol=p^2)
  E_i <- matrix(1:p^2,ncol=p)
  e <- E_i[lower.tri(E_i)]
  for(i in 1:nrow(E_l)){
    E_l[i,][e[i]]=1
  }
  return(E_l)
}

up.elim <- function(p){
  E_u <- matrix(0,nrow=p*(p-1)/2,ncol=p^2)
  E_i <- matrix(1:p^2,ncol=p)
  e <- E_i[upper.tri(E_i)]
  for(i in 1:nrow(E_u)){
    E_u[i,][e[i]]=1
  }
  return(E_u)
}

diag.elim <- function(p){
  E_d <- matrix(0,nrow=p,ncol=p^2)
  E_i <- matrix(1:p^2,ncol=p)
  e <- E_i[diag(E_i)]
  for(i in 1:nrow(E_d)){
    E_d[i,][e[i]]=1
  }
  return(E_d)
}

lo.dupl <- function(p){
  E_p <- matrix(0,nrow=p^2,ncol= p*(p-1)/2)
  E_i <- matrix(1:p^2,ncol=p)
  e <- E_i[lower.tri(E_i)]
  for(i in e){
    E_p[i,][which(e==i)]=1
  }
  return(E_p)
}

### dvecK/dvec(logK)
dvecl <- function(B,s=svd(B)){
  #s <- svd(B)
  p <- nrow(B)
  l <- s$d
  Q <- s$u
  P <- matrix(0,ncol=p^2,nrow=p^2)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j){
        P[((i-1)*p+j),((i-1)*p+j)]= exp(l[i])
      } else {
        P[((i-1)*p+j),((i-1)*p+j)] = (exp(l[i])-exp(l[j]))/(l[i]-l[j])
      }
    }
  }
  P1 <- kronecker(Q,Q) %*% P %*% t(kronecker(Q,Q))
  return(P1)
}

# inverse of dvecl
inv.dvecl <- function(B,s=svd(B)){
  #s <- svd(B)
  p <- nrow(B)
  l <- s$d
  Q <- s$u
  Qinv <- t(Q)
  P <- matrix(0,ncol=p^2,nrow=p^2)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j){
        P[((i-1)*p+j),((i-1)*p+j)]= exp(l[i])
      } else {
        P[((i-1)*p+j),((i-1)*p+j)] = (exp(l[i])-exp(l[j]))/(l[i]-l[j])
      }
    }
  }
  P1 <- kronecker(Qinv,Qinv) %*% diag(1/diag(P)) %*% t(kronecker(Qinv,Qinv))
  return(P1)
}


beta.elim <- function(p){
  E_b <- matrix(0,nrow=(p-1)*p, ncol=p^2)
  E_i <- matrix(1:p^2,ncol=p)
  e <- c(E_i[-1,1],as.numeric(E_i[-1,-1]))
  for(i in 1:nrow(E_b)){
    E_b[i,][e[i]]= 1
  }
  return(E_b)
}

log_m <- function(B,s=svd(B)){
  p <- nrow(B)
  l <- s$d
  Q <- s$u
  return(Q %*% diag(log(l))  %*% t(Q))
}

# NOTE: named deriv_b (not deriv.b) so roxygen2 does not mistake it for an
# S3 method of the stats::deriv generic. Internal helper used by sgclm().
deriv_b <- function(A,dep.v =1){
  Rinv <- as.matrix(solve(A[-dep.v,-dep.v]))
  p <- ncol(A)
  deriv.beta1 <- Rinv
  deriv.beta2 <- -(kronecker(Rinv,Rinv)) %*% (kronecker(diag(p-1),A[dep.v,-dep.v]))
  # deriv.beta1 <- solve(Rinv)
  # deriv.beta2 <- - kronecker(t(Rinv %*% R1[-dep.v,dep.v]), Rinv)
  deriv.beta <- rbind(deriv.beta1,deriv.beta2)
  return(t(deriv.beta))
}


