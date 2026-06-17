# dev_smoke_test.R --------------------------------------------------------
# Manual smoke tests for the cleaned-up SGCTools dev package.
# Run after: devtools::document(); devtools::load_all()
# Not part of the package build (see .Rbuildignore).
# -------------------------------------------------------------------------

set.seed(1)
n <- 500; m <- 15
tp <- seq(0, 1, length = m)
d  <- abs(outer(tp, tp, "-"))
Sigma <- exp(-2 * d)                      # latent correlation surface
y <- mvtnorm::rmvnorm(n, sigma = Sigma)   # latent Gaussian process

# helper: report quietly, check covariance dims + leading eigenvalues
report <- function(label, ff) {
  stopifnot(is.matrix(as.matrix(ff$cov)),
            nrow(as.matrix(ff$cov)) == m,
            length(ff$evalues) == 4,
            all(ff$evalues >= -1e-6))     # PD up to tolerance
  cat(sprintf("%-26s OK | evalues: %s\n", label,
              paste(round(ff$evalues, 3), collapse = ", ")))
  invisible(ff)
}

cat("=== fpca.sgc.lat across all four types ===\n")

## 1. continuous -----------------------------------------------------------
z_cont <- y
report("cont", fpca.sgc.lat(X = z_cont, type = "cont", argvals = tp, df = 4))

## 2. binary ---------------------------------------------------------------
z_bin <- (y > 0.5) * 1
report("bin", fpca.sgc.lat(X = z_bin, type = "bin", argvals = tp, df = 4))

## 3. ordinal (4 categories: 0..3, common cutpoints across all columns) ----
brks <- c(-Inf, -0.6, 0.1, 0.6, Inf)
z_ord <- apply(y, 2, function(col) as.numeric(cut(col, breaks = brks)) - 1)
# confirm every column actually shows all 4 categories (ord branch assumes this)
stopifnot(all(apply(z_ord, 2, function(c) length(unique(c))) == 4))
report("ord", fpca.sgc.lat(X = z_ord, type = "ord", argvals = tp, df = 4))

## 4. truncated (zero-inflated continuous) ---------------------------------
z_trunc <- ifelse(y > 0, y, 0)            # point mass at 0 + positive part
cat(sprintf("trunc zero-fraction: %.2f\n", mean(z_trunc == 0)))
report("trunc", fpca.sgc.lat(X = z_trunc, type = "trunc", argvals = tp, df = 4))

## 5. NEW input guard: invalid type should error cleanly -------------------
cat("\n=== type validation guard ===\n")
res <- tryCatch(fpca.sgc.lat(X = z_bin, type = "binary"),
                error = function(e) conditionMessage(e))
cat("invalid type ->", res, "\n")
stopifnot(grepl("must be one of", res))

## 6. scores + missing data path (exercises getLatentPreds / recover_row) --
cat("\n=== scores with missing data (impute) ===\n")
z_miss <- z_bin
z_miss[sample(length(z_miss), 0.15 * length(z_miss))] <- NA   # 15% MCAR
ff_sc <- fpca.sgc.lat(X = z_miss, type = "bin", argvals = tp, df = 4,
                      scores = TRUE, impute = TRUE)
stopifnot(nrow(ff_sc$scores) == n, ncol(ff_sc$scores) == 4)
cat("scores dim:", paste(dim(ff_sc$scores), collapse = " x "), "OK\n")

# -------------------------------------------------------------------------
cat("\n=== scalar mixed-data pipeline (fromXtoRMixed / sgclm / getLatentPreds) ===\n")

p <- 5
S <- clusterGeneration::rcorrmatrix(p)
L <- mvtnorm::rmvnorm(200, sigma = S)
# C, B, O, O, T
X <- data.frame(
  x1 = L[, 1]^3,
  y1 = as.numeric(L[, 2] >= 1),
  x2 = as.numeric(cut(L[, 3], breaks = c(-Inf, -0.2, 0.1, 0.3, Inf))) - 1,
  x3 = as.numeric(cut(L[, 4], breaks = c(-Inf, -0.3, 0.05, 0.2, Inf))) - 1,
  x4 = ifelse(L[, 5] > 0, L[, 5], 0)
)

K <- Kendall_mixed(X)
stopifnot(all(abs(diag(K) - 1) < 1e-12))
cat("Kendall_mixed: diag == 1  OK\n")

R <- fromXtoRMixed(X, use.nearPD = TRUE)
stopifnot(isSymmetric(unname(R$hatR)), all(abs(diag(R$hatR) - 1) < 1e-6))
cat("fromXtoRMixed: symmetric corr matrix, unit diag  OK\n")

pc <- princomp(covmat = R$hatR, cor = TRUE)
cat("princomp on hatR: top sdev =", round(pc$sdev[1], 3), "OK\n")

lm1 <- sgclm(y1 ~ x1 + x2 + x4, data = X)
cat("sgclm coefficient table:\n"); print(lm1$coef)

Lhat <- getLatentPreds(X, lat.cov = R$hatR)
stopifnot(nrow(Lhat) == nrow(X), ncol(Lhat) == ncol(X))
cat("getLatentPreds dim:", paste(dim(Lhat), collapse = " x "), "OK\n")

cat("\nAll smoke tests passed.\n")
