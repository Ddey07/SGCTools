require(SGCTools)

matern <- function (u, phi, kappa)
{
  if (is.vector(u))
    names(u) <- NULL
  if (is.matrix(u))
    dimnames(u) <- list(NULL, NULL)
  uphi <- u/phi
  uphi <- ifelse(u > 0,
                 ((2^(-(kappa - 1))) / ifelse(0, Inf, gamma(kappa))) *
                   (uphi^kappa) * besselK(x = uphi, nu = kappa),
                 1)
  uphi[u > 600 * phi] <- 0
  return(uphi)
}

n <- 1000
m <- 15
delta <- 0.5 # cutoff
cmb <- combn(m,2)
tp = seq(0,1,length=m)
d = abs(outer(tp,tp,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
phi=2 # length scale
l=0.01

# Generate covariance matrix from stationary kernel
Sigma_SE = matern(d,phi= 1/phi, kappa= 3.5) # Matern exponential kernel

# generate latent process
set.seed(1)
y = mvtnorm::rmvnorm(n,sigma=Sigma_SE)

# generate observed process based on the cutoff
z = y
z[z>delta]=1
z[z<=delta]=0

# run fpca.sgc.lat to get functional principal component analysis of binary data
ff = fpca.sgc.lat(X=z,type="bin",argvals = tp, df= 4)

# compare truth vs estimate
par(mfrow=c(1,2))
image(as.matrix(Sigma_SE), main="Truth")
image(as.matrix(ff$cov), main="Estimate")
