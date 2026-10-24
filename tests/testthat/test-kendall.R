tau_a_brute <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; n <- length(x); s <- 0
  for (i in 1:(n - 1)) for (j in (i + 1):n) s <- s + sign((x[i] - x[j]) * (y[i] - y[j]))
  s / choose(n, 2)
}

test_that("Kendall_mixed returns 1 for perfectly concordant data", {
  for (n in c(2, 5, 10, 11, 20, 37)) {
    x <- seq_len(n)
    expect_equal(Kendall_mixed(cbind(x, 2 * x + 1))[1, 2], 1)
    expect_equal(Kendall_mixed(cbind(x, -x))[1, 2], -1)
  }
})

test_that("Kendall_mixed equals tau_a with ties and missing values", {
  set.seed(1)
  for (r in 1:60) {
    n <- sample(3:80, 1)
    z <- matrix(rnorm(2 * n), n) %*% chol(matrix(c(1, 0.6, 0.6, 1), 2))
    a <- switch(r %% 4 + 1, z[, 1], (z[, 1] > 0) * 1, pmax(z[, 1], 0), findInterval(z[, 1], c(-0.5, 0.5)))
    b <- switch((r %/% 4) %% 4 + 1, z[, 2], (z[, 2] > 0.3) * 1, pmax(z[, 2], 0), findInterval(z[, 2], c(-1, 0)))
    a[sample(n, floor(0.2 * n))] <- NA; b[sample(n, floor(0.2 * n))] <- NA
    if (sum(!is.na(a) & !is.na(b)) < 2) next
    expect_equal(Kendall_mixed(cbind(a, b))[1, 2], tau_a_brute(a, b), tolerance = 1e-12)
  }
})

test_that("bridgeF_tt and fromKtoR_mixed keep named arguments for Vectorize", {
  expect_identical(names(formals(SGCTools:::bridgeF_tt)), c("r", "zratio1", "zratio2"))
  v <- Vectorize(SGCTools:::bridgeF_tt, vectorize.args = c("r", "zratio1", "zratio2"))
  expect_length(v(c(0.2, 0.5), c(0.3, 0.3), c(0.4, 0.4)), 2)
  expect_true(all(c("K12", "zratio1", "zratio2", "type1", "type2", "tol") %in%
                  names(formals(SGCTools:::fromKtoR_mixed))))
})
