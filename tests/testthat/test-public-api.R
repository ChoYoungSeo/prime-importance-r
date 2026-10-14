make_sparse_example <- function(n = 300L) {
  set.seed(5)
  X <- matrix(runif(n * 8L, -1, 1), nrow = n)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  y <- 2 * X[, 1] - X[, 2] + 1.5 * X[, 4]^2 +
    rnorm(n, sd = 0.1)
  list(X = X, y = y)
}

test_that("locked defaults match the manuscript specification", {
  constants <- prime:::.prime_constants
  expect_identical(constants$ridge_cap_multiplier, 1.5)
  expect_identical(constants$neighborhood_a, 1.5)
  expect_identical(constants$neighborhood_b, 3)
  expect_identical(constants$cap_fraction, 0.5)
  expect_identical(constants$ridge_lambda, 1e-6)
  expect_identical(prime:::.prime_max_mask_dimension, 30L)
  expect_false("alpha" %in% names(formals(prime)))
  expect_false("alpha" %in% names(formals(prime_screen)))

  for (n in c(100L, 500L, 2000L)) {
    for (d in c(1L, 3L, 5L)) {
      expected <- ceiling(max(1.5 * sqrt(n), 3 * d))
      cap <- min(n - 1L, floor(n / 2L))
      if (expected <= d + 1L || expected > cap) expected <- NA_integer_
      expect_identical(
        prime:::.prime_neighborhood_size(n, d),
        as.integer(expected)
      )
    }
  }
  expect_identical(prime:::.prime_dimension_limit(100L, 100L), 16L)
  old_threads <- Sys.getenv("PRIME_NUM_THREADS", unset = NA_character_)
  on.exit({
    if (is.na(old_threads)) {
      Sys.unsetenv("PRIME_NUM_THREADS")
    } else {
      Sys.setenv(PRIME_NUM_THREADS = old_threads)
    }
  }, add = TRUE)
  Sys.setenv(PRIME_NUM_THREADS = "8")
  expect_identical(prime:::.prime_worker_count(10L), 8L)
})

test_that("screen and PRIME allocation expose a coherent public result", {
  example <- make_sparse_example()
  fit <- prime(example$X, example$y)

  expect_s3_class(fit, "prime_fit")
  expect_true(all(fit$selected %in% fit$screen$candidates))
  expect_identical(
    names(fit$effects),
    colnames(example$X)
  )
  expect_true(all(c("risk_full", "noise_floor", "signal") %in% names(fit)))
  expect_true(all(c(
    "complement_risk", "adjusted_risk"
  ) %in% names(fit$subset_values)))
  expect_equal(fit$risk_full, fit$noise_floor, tolerance = 1e-14)
  expect_equal(
    fit$subset_values$adjusted_risk,
    pmax(fit$subset_values$complement_risk - fit$noise_floor, 0),
    tolerance = 1e-12
  )
  expect_false("tau_full" %in% names(fit))
  expect_false("complement_tau" %in% names(fit$subset_values))
  expect_equal(sum(fit$effects), 1, tolerance = 1e-8)
  expect_true(all(fit$effects[-fit$selected] == 0))
  expect_false("alpha" %in% names(fit$settings))
  expect_false("backward_cutoff" %in% names(fit$settings))
  expect_identical(fit$settings$lambda, 1e-6)
  expect_identical(fit$settings$local_lambda, 1e-6)
  expect_identical(fit$settings$ridge_lambda, 1e-6)
  expect_identical(fit$settings$residual_df, "K-d-1")
})

test_that("fit and screen summaries are structured and printable", {
  example <- make_sparse_example()
  fit <- prime(example$X, example$y)

  fit_summary <- summary(fit)
  screen_summary <- summary(fit$screen)

  expect_s3_class(fit_summary, "summary_prime_fit")
  expect_identical(fit_summary$n, nrow(example$X))
  expect_identical(fit_summary$p, ncol(example$X))
  expect_equal(sum(fit_summary$importance[, "Estimate"]), 1, tolerance = 1e-8)
  expect_true(is.call(fit_summary$call))
  expect_true(is.finite(fit_summary$explained_fraction))
  expect_s3_class(screen_summary, "summary_prime_screen")
  expect_identical(screen_summary$retained, fit$selected_names)
  expect_identical(
    screen_summary$backward_evaluations,
    nrow(fit$screen$backward$trace)
  )
  expect_identical(
    screen_summary$backward_deletions,
    sum(fit$screen$backward$trace$removed)
  )
  expect_output(print(fit), "PRIME fit")
  expect_output(print(fit_summary), "PRIME importance")
  expect_output(print(fit$screen), "PRIME screen")
  expect_output(print(screen_summary), "Forward selection")
})

test_that("plot method accepts a fitted PRIME object", {
  example <- make_sparse_example()
  fit <- prime(example$X, example$y)
  output <- tempfile(fileext = ".png")

  grDevices::png(output)
  plotted <- plot(fit)
  grDevices::dev.off()

  expect_identical(plotted, fit)
  expect_gt(file.info(output)$size, 0)
})

test_that("local slope nugget does not tune the global ridge seed", {
  example <- make_sparse_example()
  default_data <- prime:::.prime_prepare(
    example$X, example$y, TRUE, NULL, 1e-6
  )
  altered_data <- prime:::.prime_prepare(
    example$X, example$y, TRUE, NULL, 1
  )
  default_ridge <- prime:::.prime_ridge_initial(default_data)
  altered_ridge <- prime:::.prime_ridge_initial(altered_data)

  expect_identical(default_ridge$selected, altered_ridge$selected)
  expect_equal(default_ridge$weights, altered_ridge$weights)
})

test_that("the calibrated path starts from the BIC ridge seed", {
  example <- make_sparse_example()
  screen <- prime_screen(example$X, example$y)
  forward <- screen$forward

  expect_lte(length(screen$candidates), forward$dimension_limit)
  expect_true(all(forward$ridge_initial %in% screen$candidates))
  expect_identical(
    length(forward$ridge_initial), forward$ridge_chosen_size
  )
  expect_identical(
    forward$ridge_chosen_size,
    forward$ridge_path$size[which.min(forward$ridge_path$bic)]
  )
  expect_setequal(
    screen$candidates,
    c(forward$ridge_initial, forward$local_path)
  )
  expect_identical(
    length(forward$local_path),
    sum(forward$trace$accepted)
  )
  expect_true(all(
    forward$trace$calibrated_gain[forward$trace$accepted] > 0
  ))
  expect_equal(
    forward$trace$threshold,
    forward$trace$query_scale *
      sqrt(2 * log(pmax(2L, forward$trace$dimension))),
    tolerance = 1e-12
  )
  expect_equal(
    forward$trace$calibrated_gain,
    pmax(forward$trace$raw_gain - forward$trace$threshold, 0),
    tolerance = 1e-12
  )
  if (forward$stop_reason ==
      "raw-gain leader's calibrated gain was nonpositive") {
    expect_false(tail(forward$trace$accepted, 1L))
    expect_equal(tail(forward$trace$calibrated_gain, 1L), 0)
  }
})

test_that("backward trace records the querywise calibration", {
  example <- make_sparse_example()
  trace <- prime_screen(example$X, example$y)$backward$trace

  expect_gt(nrow(trace), 0L)
  expect_equal(
    trace$threshold,
    trace$query_scale * sqrt(2 * log(pmax(2L, trace$dimension))),
    tolerance = 1e-12
  )
  expect_equal(trace$clipped_loss, pmax(trace$raw_loss, 0))
  expect_equal(
    trace$calibrated_loss,
    pmax(trace$raw_loss - trace$threshold, 0),
    tolerance = 1e-12
  )
  expect_identical(trace$removed, trace$calibrated_loss == 0)
})

test_that("exact proportional allocation is not capped at 12 coordinates", {
  dimension <- 13L
  game <- seq(0, 1, length.out = bitwShiftL(1L, dimension))
  effects <- prime_allocate(game)

  expect_length(effects, dimension)
  expect_equal(sum(effects), 1, tolerance = 1e-10)
})

test_that("correlated suppressor variables survive ridge-start screening", {
  set.seed(5)
  n <- 500L
  p <- 10L
  X <- matrix(0, nrow = n, ncol = p)
  X[, 1] <- runif(n, -1, 1)
  X[, 2] <- runif(n, -1, 1)
  X[, 3] <- 0.9 * X[, 2] + 0.1 * rnorm(n)
  X[, 4] <- 0.9 * X[, 1] + 0.1 * rnorm(n)
  X[, 5:p] <- matrix(runif(n * (p - 4L), -1, 1), nrow = n)
  y <- X[, 1] + X[, 2] - X[, 3] + X[, 5]^2 +
    0.1 * rnorm(n)

  screen <- prime_screen(X, y)
  expect_true(all(c(1L, 2L, 3L, 5L) %in% screen$candidates))
  expect_true(all(c(1L, 2L, 3L, 5L) %in% screen$selected))
})

test_that("ridge-start screening respects local-fit feasibility", {
  set.seed(204709101L)
  n <- 100L
  p <- 100L
  X <- stats::pnorm(MASS::mvrnorm(
    n, mu = numeric(p), Sigma = diag(p)
  ))
  signal <- X[, 1] + X[, 2] + 3 * X[, 1] * X[, 3]
  y <- signal + stats::rnorm(n, sd = 0.5)

  screen <- prime_screen(X, y)
  expect_lte(length(screen$candidates), screen$forward$dimension_limit)
})

test_that("invalid input fails clearly", {
  example <- make_sparse_example()
  bad <- example$X
  bad[1, 1] <- NA_real_
  expect_error(prime_screen(bad, example$y), "finite")
  expect_error(
    prime_screen(example$X, rep(1, nrow(example$X))),
    "nonconstant"
  )
  expect_error(
    prime_screen(example$X, example$y, twin_r = nrow(example$X)),
    "floor\\(n/2\\)"
  )
  constant <- cbind(example$X, constant = 1)
  expect_error(
    prime_screen(constant, example$y),
    "must not contain constant columns",
    fixed = TRUE
  )
})
