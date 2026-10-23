## Tests for the getLatentPreds() / preprocess_data() fix (SGCTools 1.1.4).
## Small problems with an exchangeable latent correlation keep these fast.

lat_cor <- function(p, rho = 0.5) { R <- matrix(rho, p, p); diag(R) <- 1; R }
pre <- function(...) SGCTools:::preprocess_data(...)

## the pre-1.1.4 cutoff construction, used to check backward compatibility
old_ord_cutoff <- function(x) {
  cats <- sort(unique(stats::na.omit(x)))
  c(-Inf, unlist(lapply(cats, function(l) stats::qnorm(1 - mean(x >= l, na.rm = TRUE))))[-1], Inf)
}

set.seed(1)
n <- 120; p <- 6
Z <- function(n, R) matrix(stats::rnorm(n * ncol(R)), n) %*% chol(R)
Zs <- Z(n, lat_cor(p))
make_ord <- function(z, cuts = c(-0.6, 0.1, 0.6)) findInterval(z, cuts)   # codes 0..3

test_that("a sparse ordinal column missing an interior level no longer gives NA bounds", {
  X <- cbind(apply(Zs[, 1:3], 2, make_ord), Zs[, 4:6])
  X[X[, 1] == 1, 1] <- NA               # level 1 absent from column 1, levels 0, 2, 3 present
  X[sample(n, 60), 2] <- NA             # extra sparsity
  type <- c(rep("ord", 3), rep("cont", 3))
  L <- getLatentPreds(X, type = type, lat.cov.est = lat_cor(p), impute.missing = TRUE)
  expect_true(all(is.finite(L)))
  cut1 <- pre(X, type = type)$features[[1]]$cutoff
  expect_length(cut1, 3 + 2)            # -Inf, three cutoffs for levels 1..3, Inf
  expect_equal(cut1[2], cut1[3])        # absent level 1 -> zero-width interval
})

test_that("results are identical to the old code when every level is present", {
  X <- cbind(apply(Zs[, 1:3], 2, make_ord), Zs[, 4:6])
  type <- c(rep("ord", 3), rep("cont", 3))
  expect_true(all(apply(X[, 1:3], 2, function(x) setequal(unique(x), 0:3))))
  for (j in 1:3) expect_equal(pre(X, type = type)$features[[j]]$cutoff, old_ord_cutoff(X[, j]))
})

test_that("the supplied type is used, not a guess from the data", {
  X <- cbind(apply(Zs[, 1:2], 2, make_ord), Zs[, 3:6])
  X[, 1] <- pmin(X[, 1], 1)             # ordinal column that only shows codes 0/1
  type <- c("ord", "ord", rep("cont", 4))
  expect_equal(pre(X, type = type)$type, type)
  expect_true(all(is.finite(getLatentPreds(X, type = type, lat.cov.est = lat_cor(p)))))
})

test_that("ordinal components with different numbers of levels can be mixed", {
  X <- cbind(make_ord(Zs[, 1]), findInterval(Zs[, 2], c(-1, -0.5, 0, 0.5, 1, 1.5)), Zs[, 3:6])  # 4 and 7 levels
  type <- c("ord", "ord", rep("cont", 4))
  f <- pre(X, type = type)$features
  expect_equal(f[[1]]$cutoff[1:4], old_ord_cutoff(X[, 1])[1:4])   # extra levels of the 4-level column are +Inf
  expect_true(all(is.infinite(f[[1]]$cutoff[5:8])))
  expect_true(all(is.finite(getLatentPreds(X, type = type, lat.cov.est = lat_cor(p)))))
})

test_that("an ordinal column with no 0 observed gets correct lower bounds", {
  X <- cbind(make_ord(Zs[, 1]), Zs[, 2:6]); X[X[, 1] == 0, 1] <- NA
  cut <- pre(X, type = c("ord", rep("cont", 5)))$features[[1]]$cutoff
  expect_equal(cut[2], -Inf)            # P(X >= 1) = 1 among observed
  expect_true(all(is.finite(getLatentPreds(X, type = c("ord", rep("cont", 5)), lat.cov.est = lat_cor(p)))))
})

test_that("degenerate columns are handled: all-NA continuous/ordinal, constant binary", {
  X <- cbind(make_ord(Zs[, 1]), as.numeric(Zs[, 2] > 0), Zs[, 3:6])
  X[, 2] <- 0                           # binary column with no 1 observed
  X[, 4] <- NA                          # continuous column never observed
  type <- c("ord", "bin", rep("cont", 4))
  f <- pre(X, type = type)$features
  expect_null(f[[4]]$fn)
  L <- getLatentPreds(X, type = type, lat.cov.est = lat_cor(p), impute.missing = TRUE)
  expect_true(all(is.finite(L)))
  X2 <- X; X2[, 1] <- NA                # ordinal column never observed
  expect_true(all(is.finite(getLatentPreds(X2, type = type, lat.cov.est = lat_cor(p)))))
})

test_that("truncated columns still work, including with no zero observed", {
  tr <- Zs[, 1] * (Zs[, 1] > 0.5)
  X <- cbind(tr, Zs[, 2:6]); type <- c("trunc", rep("cont", 5))
  expect_true(all(is.finite(getLatentPreds(X, type = type, lat.cov.est = lat_cor(p)))))
  X[X[, 1] == 0, 1] <- NA
  expect_true(all(is.finite(getLatentPreds(X, type = type, lat.cov.est = lat_cor(p)))))
})

test_that("type = NULL keeps the automatic type detection", {
  X <- cbind(make_ord(Zs[, 1]), as.numeric(Zs[, 2] > 0), Zs[, 3:6])
  expect_equal(pre(X)$type, c("ord", "bin", rep("cont", 4)))
  expect_true(all(is.finite(getLatentPreds(X, lat.cov.est = lat_cor(p)))))
})

test_that("invalid inputs give informative errors", {
  X <- cbind(make_ord(Zs[, 1]), Zs[, 2:6])
  expect_error(getLatentPreds(X, type = c("ord", "cont"), lat.cov.est = lat_cor(p)), "length ncol")
  expect_error(getLatentPreds(X, type = c("ordinal", rep("cont", 5)), lat.cov.est = lat_cor(p)), "unknown 'type'")
  Xn <- X; Xn[1, 1] <- -1
  expect_error(getLatentPreds(Xn, type = c("ord", rep("cont", 5)), lat.cov.est = lat_cor(p)), "non-negative integers")
  Xf <- X; Xf[1, 1] <- 1.5
  expect_error(getLatentPreds(Xf, type = c("ord", rep("cont", 5)), lat.cov.est = lat_cor(p)), "non-negative integers")
  Xb <- cbind(as.numeric(Zs[, 1] > 0), Zs[, 2:6]); Xb[1, 1] <- 2
  expect_error(getLatentPreds(Xb, type = c("bin", rep("cont", 5)), lat.cov.est = lat_cor(p)), "0/1")
})

test_that("fpca.sgc.lat latent scores run on sparse ordinal functional data", {
  set.seed(2); m <- 8; nn <- 150
  tt <- seq(0, 1, length.out = m)
  C <- exp(-abs(outer(tt, tt, "-")) / 0.4)
  V <- matrix(stats::rnorm(nn * m), nn) %*% chol(C)
  X <- apply(V, 2, make_ord)
  X[matrix(stats::runif(nn * m) < 0.6, nn)] <- NA       # 40% observed
  X[1, ] <- NA; X[1, 3] <- make_ord(V[1, 3])            # a subject with exactly one observation
  X[2, ] <- NA                                          # a subject with none
  fit <- suppressWarnings(fpca.sgc.lat(X, type = "ord", df = 4, npc = 2, scores = TRUE, min_no_pairs = 10))
  obs <- !is.na(X)
  expect_true(all(is.finite(fit$latent[obs])))          # impute = FALSE: predictions at observed cells
  has_obs <- rowSums(obs) > 0
  expect_true(all(is.finite(fit$scores[has_obs, ])))
  expect_true(all(is.na(fit$scores[!has_obs, ])))       # no data -> NA scores, not an error
})

test_that("data.frame input still works (fromXtoRMixed path)", {
  df <- data.frame(a = Zs[, 1], b = Zs[, 2] > 0, c = make_ord(Zs[, 3]))
  expect_silent(pre(df))
  expect_equal(pre(df)$type, c("cont", "bin", "ord"))
})

test_that("fpca.sgc.lat (ordinal) handles columns that miss a level, including column 1", {
  set.seed(4); m <- 6; nn <- 300
  tt <- seq(0, 1, length.out = m)
  V <- matrix(stats::rnorm(nn * m), nn) %*% chol(exp(-abs(outer(tt, tt, "-")) / 0.5))
  X <- apply(V, 2, make_ord)
  X[X[, 1] == 3, 1] <- NA               # column 1 misses the top level (used to fix the level count)
  X[X[, 4] == 1, 4] <- NA               # column 4 misses an interior level
  fit <- suppressWarnings(fpca.sgc.lat(X, type = "ord", df = 4, npc = 2, min_no_pairs = 10))
  expect_true(all(is.finite(as.matrix(fit$cov))))
  expect_equal(dim(as.matrix(fit$cov)), c(m, m))
})

test_that("fpca.sgc.lat rejects badly coded ordinal data", {
  X <- matrix(c(1.5, rep(0:3, 30)), ncol = 1)[1:120, , drop = FALSE]
  X <- cbind(X, X, X)
  expect_error(fpca.sgc.lat(X, type = "ord", df = 4), "non-negative integers")
})
