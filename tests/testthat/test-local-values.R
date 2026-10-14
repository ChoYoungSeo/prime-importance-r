make_local_data <- function(n = 120L, p = 8L) {
  set.seed(20260816)
  X <- matrix(rnorm(n * p), nrow = n)
  y <- 2 * X[, 1] - X[, 3]^2 + sin(X[, 5]) + rnorm(n, sd = 0.2)
  prime:::.prime_prepare(
    X, y,
    rescale = TRUE,
    twin_r = NULL,
    lambda = 1e-6
  )
}

test_that("compiled querywise risks equal direct local residualization", {
  data <- make_local_data()
  subset <- c(1L, 3L, 5L)
  k <- prime:::.prime_neighborhood_size(data$n, length(subset))
  compiled <- prime:::.prime_variance_details(data, subset, k)
  direct <- prime:::.prime_variance_details_r(data, subset, k)

  expect_equal(compiled$mean_sse, direct$mean_sse, tolerance = 1e-10)
  expect_equal(compiled$residual_df, direct$residual_df)
  expect_equal(compiled$risk, direct$risk, tolerance = 1e-10)
  expect_equal(
    compiled$contribution, direct$contribution, tolerance = 1e-10
  )
  expect_equal(mean(compiled$contribution), compiled$risk, tolerance = 1e-12)
})

test_that("batched forward scores reproduce raw risk rankings", {
  data <- make_local_data(p = 12L)
  selected <- c(1L, 3L, 5L)
  candidates <- setdiff(seq_len(data$p), selected)
  k <- prime:::.prime_neighborhood_size(
    data$n, length(selected) + 1L
  )
  expected <- vapply(candidates, function(candidate) {
    prime:::.prime_forward_score(data, c(selected, candidate), k)
  }, numeric(1))
  risks <- vapply(candidates, function(candidate) {
    prime:::.prime_risk(data, c(selected, candidate), k)
  }, numeric(1))
  actual <- prime:::.prime_forward_values(
    data, selected, candidates, k
  )

  expect_equal(actual, expected, tolerance = 1e-11)
  expect_identical(which.max(actual), which.min(risks))
})

test_that("forward ranking is preserved below the zero-signal boundary", {
  set.seed(85)
  X <- matrix(rnorm(15L * 8L), nrow = 15L)
  y <- rnorm(15L)
  data <- prime:::.prime_prepare(X, y, TRUE, NULL, 1e-6)
  selected <- 1L
  candidates <- 2:8
  k <- prime:::.prime_neighborhood_size(data$n, 2L)
  scores <- prime:::.prime_forward_values(data, selected, candidates, k)
  risks <- vapply(candidates, function(candidate) {
    prime:::.prime_risk(data, c(selected, candidate), k)
  }, numeric(1))

  expect_true(all(scores < 0))
  expect_identical(which.max(scores), which.min(risks))
})

test_that("universal deletion diagnostic uses paired querywise losses", {
  data <- make_local_data(n = 180L)
  selected <- c(1L, 3L, 5L)
  variable <- 5L
  k <- prime:::.prime_neighborhood_size(data$n, length(selected))
  full <- prime:::.prime_variance_details(data, selected, k)
  reduced <- prime:::.prime_variance_details(
    data, setdiff(selected, variable), k
  )
  paired <- reduced$contribution - full$contribution
  actual <- prime:::.prime_deletion_diagnostic(
    data, selected, variable, k, full_details = full
  )

  expect_equal(actual$paired_loss, paired, tolerance = 1e-12)
  expect_equal(actual$raw_loss, mean(paired), tolerance = 1e-12)
  expect_equal(
    actual$query_scale,
    stats::sd(paired) / sqrt(length(paired)),
    tolerance = 1e-12
  )
  expect_equal(
    actual$threshold,
    actual$query_scale * sqrt(2 * log(length(selected))),
    tolerance = 1e-12
  )
  expect_equal(
    actual$calibrated_loss,
    max(actual$raw_loss - actual$threshold, 0),
    tolerance = 1e-12
  )
})

test_that("forward diagnostic uses paired querywise risk gains", {
  data <- make_local_data(n = 180L)
  selected <- c(1L, 3L)
  variable <- 5L
  dimension <- length(selected) + 1L
  k <- prime:::.prime_neighborhood_size(data$n, dimension)
  current <- prime:::.prime_variance_details(data, selected, k)
  augmented <- prime:::.prime_variance_details(
    data, c(selected, variable), k
  )
  paired <- current$contribution - augmented$contribution
  actual <- prime:::.prime_forward_diagnostic(
    data, selected, variable, k
  )

  expect_equal(actual$paired_gain, paired, tolerance = 1e-12)
  expect_equal(actual$raw_gain, mean(paired), tolerance = 1e-12)
  expect_equal(
    actual$threshold,
    actual$query_scale * sqrt(2 * log(dimension)),
    tolerance = 1e-12
  )
  expect_equal(
    actual$calibrated_gain,
    max(actual$raw_gain - actual$threshold, 0),
    tolerance = 1e-12
  )
})

test_that("all-query nearest neighbors agree with FNN", {
  data <- make_local_data()
  subset <- c(1L, 2L, 4L)
  Xs <- data$X[, subset, drop = FALSE]
  k <- 20L
  actual <- prime:::.prime_query_neighbors(
    Xs, seq_len(nrow(Xs)), k
  )
  expected <- FNN::get.knn(Xs, k = k, algorithm = "kd_tree")$nn.index
  expect_identical(actual, expected)
})
