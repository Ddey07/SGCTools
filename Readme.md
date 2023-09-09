<!-- README.md is generated from README.Rmd. Please edit that file -->

# SGCTools: Estimation tools for multivariate mixed data and continous/truncated/discrete functional data using Semiparametric Gaussian Copula (SGC)

The R package `SGCTools` (a) estimates latent correlation for
multivariate mixed data types; continuous, truncated, ordinal (any
number of categories), and binary, (b) performs covariance estimation
and functional principal components analysis for continous, truncated,
discrete (ordinal/binary) functional data. The methods are described in
Dey, Zipunnikov (2022) <https://arxiv.org/abs/2205.06868> and Dey,
Ghosal, Merikangas, and Zipunnikov (2023)
<https://arxiv.org/abs/2306.15084>

## Installation

``` install
devtools::install_github("DDey07/SGCTools")
```

## Example

``` r
library(SGCTools)

matern <- function (u, phi, kappa)
{
  if (is.vector(u))
    names(u) <- NULL
  if (is.matrix(u))
    dimnames(u) <- list(NULL, NULL)
  uphi <- u/phi
  uphi <- ifelse(u > 0, (((2^(-(kappa - 1)))/ifelse(0, Inf,
                                                    gamma(kappa))) * (uphi^kappa) * besselK(x = uphi, nu = kappa)),
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
```

![](README-example-1.png)<!-- -->
