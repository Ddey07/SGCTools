#' @title Perform latent regression using SGC
#'
#' @description Estimate of latent regression coefficient based on outcome and covariates of mixed types (continuous/truncated/ordinal/binary) based on the Semiparametric Gaussian copula model.
#'
#' @param formula A regression formula specifying the name of the outcome and covariates in the data frame
#' @param type A vector of length p containing type of variables in \code{X}, must be in the format "cont", "trunc", "ord", "bin". If not supplied, our function tries to guess the data-type from the observations
#' @param data A numeric data frame (n by p), n is the sample size and p is the number of variables. Binary variables must be coded as 0, 1. Ordinal variables with l ordered categories must be coded as 0, 1, 2, ... , (l-1).
#' @return \code{sgclm} returns a data frame containing the latent regression coefficients, asymptotic standard error, p-values and significance indicator.
#'
#' @references
#' Dey ., Zipunnikov V. (2022) "Semiparametric Gaussian Copula Regression modeling for Mixed Data Types (SGCRM)" <https://arxiv.org/abs/2205.06868>
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
sgclm=function(formula,type=NULL,data){
  # Check if data is provided
  if (missing(data)) {
    stop("Data is missing. Please provide a data frame.")
  }

  # Check if formula is provided
  if (missing(formula)) {
    stop("Formula is missing. Please provide a formula.")
  }

  # Check if formula is valid
  if (!is.null(data) && !all(attr(terms(formula), "term.labels") %in% names(data))) {
    stop("Invalid formula. Please check the variable names in the formula.")
  }

  # Check if there are complete cases
  if (sum(complete.cases(data)) == 0) {
    stop("No complete cases in the data.")
  }

  # extract the formula components from the data
  Xd <- model.frame(formula, data = data)
  # take only complete observations
  Xd <- Xd[complete.cases(Xd),]

  # Check if there are enough observations
  if (nrow(Xd) < 2) {
    stop("Insufficient data. There must be at least two observations.")
  }

# estimate latent correlation
RR = fromXtoRMixed(Xd,deriv=TRUE,type=type)
R = RR$hatR
Rp <- RR$hatRprime
beta_est = as.numeric( solve(R[-1,-1]) %*% R[1,-1])

# dimensions of the data
p <- ncol(Xd)
n <- nrow(Xd)

# calculate variance of Kendall's Tau
CDN <- SGCTools:::conc_ties(Xd)
C <- CDN$C
D <-  CDN$D
N <- CDN$N
C_sum <- apply(C, c(2,3), function(x){sum(x, na.rm = TRUE)})
D_sum <- apply(D, c(2,3), function(x){sum(x, na.rm = TRUE)})
K1 <- (C_sum - D_sum)/(2*choose(N,2))
diag(K1) <- 1
ijs  <- cmbn <- combn(p,2)
vv = apply(ijs, 2, function(k) (C[,k[1],k[2]] - D[,k[1],k[2]])/(N[k[1],k[2]]-1) )
tauv = 4* cov(vv,use="pairwise.complete.obs")

# use delta method to calculate variance of beta
dep.v<-1
V_k <- tauv
tauv1 <-  lo.dupl(p) %*% tauv %*% t(lo.dupl(p))
svd_K = svd(K1)
logK = log_m(K1,s=svd_K)
svd_logK = svd_K
svd_logK$d = log(svd_K$d)
jacob_k = dvecl(logK, s = svd_logK)
jacob_k_inv = inv.dvecl(logK, s=svd_logK)
V_g <- lo.elim(p) %*% jacob_k_inv %*% tauv1 %*% jacob_k_inv %*% t(lo.elim(p))
diag(Rp) <- 0
delg <- lo.dupl(p) %*% diag(Rp[lower.tri(Rp)]) %*% t(lo.dupl(p))
veclk <- lo.elim(p) %*% (diag(p^2) - jacob_k %*% t(diag.elim(p)) %*% solve(diag.elim(p) %*% jacob_k %*%
                                                                             t(diag.elim(p))) %*% diag.elim(p)) %*% jacob_k %*% t(lo.elim(p) + up.elim(p))
V_s <- (delg %*% lo.dupl(p) %*% veclk) %*% V_g %*% t(delg %*% lo.dupl(p) %*% veclk)
V_sl <- lo.elim(p) %*% V_s %*% t(lo.elim(p))
deriv.beta <- deriv_b(R)
V_B <- (deriv.beta %*% beta.elim(p) %*% lo.dupl(p) %*% diag(Rp[lower.tri(Rp)]) %*% veclk) %*% V_g %*% t((deriv.beta %*% beta.elim(p) %*% lo.dupl(p) %*% diag(Rp[lower.tri(Rp)]) %*% veclk))

# create results data frame
sd_beta = sqrt(diag(V_B)/n)
tv = beta_est/sd_beta
pv = 2*(pt(abs(tv),lower.tail=FALSE,df=n-1))
signif = symnum(pv,
                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                symbols = c("***", "**", "*", ".", " "))

if(is.null(type)){
type = preprocess_data(Xd)$type}
res_df = data.frame("Type" = type[-1], "Estimate" = beta_est, "Std.Error" = sd_beta, 't value' = tv, "Pr(>|t|)" = pv, "Signif"= signif)

names(res_df) = c("Type", "Estimate", "Std.Error","t value","Pr(>|t|)","")
rownames(res_df) = names(Xd)[-1]

return(list(coef= res_df, legend = attr(signif,"legend")))
}


