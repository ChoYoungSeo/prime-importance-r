test_that("allocation retains positive effects below the subset zero tolerance", {
  for (tol in c(sqrt(.Machine$double.eps), 0.1)) {
    a <- tol + tol^2 / 2
    game <- c(0, a, 1, 1)
    expected <- c(a, 1) / (1 + a)

    effects <- prime_allocate(game, zero_tolerance = tol)

    expect_gt(effects[[1L]], 0)
    expect_lt(effects[[1L]], tol)
    expect_equal(as.numeric(effects), expected, tolerance = 1e-14)
    expect_equal(sum(effects), 1, tolerance = 1e-14)
    expect_equal(attr(effects, "projected_values"), game, tolerance = 1e-14)

    null_effects <- prime_allocate(c(0, tol / 2, 1, 1), zero_tolerance = tol)
    expect_equal(as.numeric(null_effects), c(0, 1), tolerance = 1e-14)
  }
})

