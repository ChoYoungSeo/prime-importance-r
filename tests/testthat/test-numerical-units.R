test_that("response units preserve a suppressor signal through ridge screening", {
  t <- rep(0:49, each = 2L)
  signal <- rep(c(-1, 1), 50L)
  X <- cbind(X1 = t + 0.1 * signal, X2 = t - 0.1 * signal)
  reference <- prime(X, signal)
  expect_identical(reference$selected, 1:2)
  expect_equal(unname(reference$effects), c(0.5, 0.5), tolerance = 1e-10)

  for (units in c(1e-100, 1e-8, -1e-8, 1e8, 1e100)) {
    fitted <- prime(X, units * signal)
    expect_identical(fitted$selected, reference$selected)
    expect_equal(fitted$effects, reference$effects, tolerance = 1e-10)
    expect_equal(fitted$screen$forward$ridge_path$bic,
                 reference$screen$forward$ridge_path$bic, tolerance = 1e-10)
    # Raw residual diagnostics, unlike the selection criterion, retain units.
    expect_equal(fitted$screen$forward$ridge_path$rss[1L] / units^2,
                 reference$screen$forward$ridge_path$rss[1L], tolerance = 1e-10)
  }
})

test_that("precomputed game allocation is homogeneous over extreme finite scales", {
  game <- c(0, 0.2, 0.8, 1)
  for (units in c(1e-250, 1e-9, 1, 1e200, 1e300)) {
    allocation <- prime_allocate(units * game)
    expect_true(all(is.finite(allocation)))
    expect_equal(as.numeric(allocation) / units, c(0.2, 0.8), tolerance = 1e-12)
    expect_equal(sum(allocation) / units, 1, tolerance = 1e-12)
    expect_equal(attr(allocation, "projected_values") / units, game,
                 tolerance = 1e-12)
  }
})

test_that("zero tolerance preserves positive games at the null boundary", {
  near_null <- c(0, 5e-11, 0, 1)
  exact <- prime_allocate(near_null, zero_tolerance = 0)
  expect_equal(as.numeric(exact), c(1, 0), tolerance = 1e-12)
  expect_identical(attr(exact, "projected_values"), near_null)
  thresholded <- prime_allocate(near_null)
  expect_equal(as.numeric(thresholded), c(0.5, 0.5), tolerance = 1e-12)

  # All proper nonempty coalitions have a tiny positive value. Symmetry and
  # efficiency determine the answer even though ordinary potentials underflow.
  tiny_game <- c(0, rep(1e-250, 6L), 1)
  tiny <- prime_allocate(tiny_game, zero_tolerance = 0)
  expect_true(all(is.finite(tiny)))
  expect_equal(as.numeric(tiny), rep(1 / 3, 3L), tolerance = 1e-12)
  expect_identical(attr(tiny, "projected_values"), tiny_game)
})

test_that("null tolerances never erase the full endpoint", {
  for (tolerance in c(1, 10)) {
    allocated <- prime_allocate(c(0, 0.2, 0.8, 1), zero_tolerance = tolerance)
    expect_equal(as.numeric(allocated), c(0.5, 0.5), tolerance = 1e-12)
    expect_identical(attr(allocated, "projected_values"), c(0, 0, 0, 1))
  }
  zero <- prime_allocate(c(0, -100, 100, 0))
  expect_identical(as.numeric(zero), c(0, 0))
  expect_identical(attr(zero, "projected_values"), rep(0, 4L))
})

test_that("invalid tolerance and unrepresentable game ratios fail clearly", {
  for (value in list(-1, Inf, NA_real_, NaN, numeric(), c(0, 0.1), "0", TRUE)) {
    expect_error(prime_allocate(c(0, 0.2, 0.8, 1), zero_tolerance = value),
                 "zero_tolerance must be a finite nonnegative numeric scalar",
                 fixed = TRUE)
  }
  expect_error(prime_allocate(c(0, 1e300, 0, 1e-300)), "numerical range")
  expect_error(prime_allocate(c(0, 1e-300, 0, 1e300)), "numerical range")
})
