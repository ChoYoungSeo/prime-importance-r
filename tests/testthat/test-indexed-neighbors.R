# Independent reference calculations deliberately avoid package neighbor helpers.
indexed_brute_neighbors <- function(X, query_index, k) {
  indices <- matrix(NA_integer_, length(query_index), k)
  distances <- matrix(NA_real_, length(query_index), k)
  ties <- logical(length(query_index))
  for (i in seq_along(query_index)) {
    query <- query_index[i] + 1L
    squared <- rowSums(sweep(X, 2L, X[query, ], "-")^2)
    eligible <- setdiff(seq_len(nrow(X)), query)
    ranked <- eligible[order(squared[eligible], eligible)]
    indices[i, ] <- ranked[seq_len(k)] - 1L
    distances[i, ] <- squared[ranked[seq_len(k)]]
    ties[i] <- k < length(ranked) &&
      squared[ranked[k]] == squared[ranked[k + 1L]]
  }
  list(indices = indices, squared_distances = distances, boundary_ties = ties)
}

indexed_brute_risk <- function(X, y, query_index, k, lambda = 1e-6) {
  neighbors <- indexed_brute_neighbors(X, query_index, k)
  sse <- vapply(seq_along(query_index), function(i) {
    index <- neighbors$indices[i, ] + 1L
    centered <- sweep(X[index, , drop = FALSE], 2L,
                      X[query_index[i] + 1L, ], "-")
    design <- cbind(1, centered)
    cross <- crossprod(design)
    diag(cross)[-1L] <- diag(cross)[-1L] + lambda
    beta <- solve(cross, crossprod(design, y[index]))
    sum((y[index] - design %*% beta)^2)
  }, numeric(1))
  df <- k - ncol(X) - 1L
  list(risk = mean(sse) / df, mean_sse = mean(sse), residual_df = df,
       contribution = sse / df, boundary_ties = neighbors$boundary_ties)
}

test_that("indexed exact neighbors match exhaustive row-index ordering", {
  set.seed(91013)
  for (dimension in c(1L, 3L, 10L)) {
    X <- matrix(rnorm(96L * dimension), nrow = 96L)
    queries <- c(95L, 0L, 31L, 0L, 58L)
    for (k in c(1L, 7L, 95L)) {
      expected <- indexed_brute_neighbors(X, queries, k)
      actual <- prime:::indexed_neighbors_cpp(X, queries, k, eps = 0)
      expect_identical(actual$indices, expected$indices)
      expect_equal(actual$squared_distances, expected$squared_distances,
                   tolerance = 1e-12)
      expect_identical(as.logical(actual$boundary_ties), expected$boundary_ties)
    }
  }
})

test_that("duplicate coordinates exclude only self and ties use row order", {
  X <- rbind(c(0, 0), c(0, 0), c(1, 0), c(-1, 0),
             c(0, 1), c(0, -1), c(2, 0))
  queries <- c(0L, 1L, 6L, 0L)
  for (k in c(1L, 2L, 3L, 6L)) {
    expected <- indexed_brute_neighbors(X, queries, k)
    actual <- prime:::indexed_neighbors_cpp(X, queries, k, eps = 0)
    expect_identical(actual$indices, expected$indices)
    expect_equal(actual$squared_distances, expected$squared_distances)
    expect_identical(as.logical(actual$boundary_ties), expected$boundary_ties)
    expect_false(any(actual$indices == queries))
  }
  identical_X <- matrix(0, nrow = 12L, ncol = 2L)
  for (eps in c(0, 0.3)) {
    actual <- prime:::indexed_neighbors_cpp(identical_X, c(0L, 7L), 5L, eps)
    expected <- indexed_brute_neighbors(identical_X, c(0L, 7L), 5L)
    expect_identical(actual$indices, expected$indices)
    expect_true(all(actual$boundary_ties))
    expect_true(all(actual$squared_distances == 0))
  }
})

test_that("ANN respects its distance guarantee and reports empirical recall", {
  set.seed(93113)
  X <- matrix(rnorm(320L * 8L), nrow = 320L)
  queries <- c(seq.int(0L, 300L, by = 17L), 0L, 319L)
  k <- 11L
  exact <- indexed_brute_neighbors(X, queries, k)
  for (eps in c(5e-8, 0.05, 0.3, 1)) {
    approximate <- prime:::indexed_neighbors_cpp(X, queries, k, eps)
    recall <- mean(vapply(seq_along(queries), function(i) {
      mean(approximate$indices[i, ] %in% exact$indices[i, ])
    }, numeric(1)))
    diagnostic <- sprintf("eps=%g, empirical recall@%d=%.3f", eps, k, recall)
    # Recall is a diagnostic: an epsilon distance bound implies no universal
    # positive lower bound on overlap with the exact top-k set.
    expect_true(is.finite(recall) && recall >= 0 && recall <= 1,
                info = diagnostic)
    expect_true(all(sqrt(approximate$squared_distances[, k]) <=
                      (1 + eps) * sqrt(exact$squared_distances[, k]) + 1e-12),
                info = diagnostic)
    expect_false(any(approximate$indices == queries), info = diagnostic)
    expect_true(all(apply(approximate$indices, 1L, anyDuplicated) == 0L),
                info = diagnostic)
    for (i in seq_along(queries)) {
      index <- approximate$indices[i, ] + 1L
      actual_squared <- rowSums(sweep(X[index, , drop = FALSE], 2L,
                                     X[queries[i] + 1L, ], "-")^2)
      expect_equal(approximate$squared_distances[i, ], actual_squared,
                   tolerance = 1e-12, info = diagnostic)
      expect_true(all(diff(approximate$squared_distances[i, ]) >= 0),
                  info = diagnostic)
    }
  }
})

test_that("streamed local risks equal exhaustive independent residualization", {
  set.seed(92113)
  X <- matrix(rnorm(73L * 3L), nrow = 73L)
  y <- 2 * X[, 1L] - X[, 2L]^2 + rnorm(nrow(X), sd = 0.3)
  queries <- c(72L, 0L, 25L, 0L, 11L, 63L, 8L, 42L, 17L)
  k <- 17L
  expected <- indexed_brute_risk(X, y, queries, k)
  for (block in c(1L, 7L, 256L)) {
    actual <- prime:::indexed_risk_cpp(
      X, y, queries, k, lambda = 1e-6, eps = 0,
      keep_contribution = TRUE, query_block_size = block
    )
    for (field in c("risk", "mean_sse", "residual_df", "contribution")) {
      expect_equal(as.numeric(actual[[field]]), as.numeric(expected[[field]]),
                   tolerance = 1e-10, info = paste("block", block, field))
    }
    expect_identical(as.logical(actual$boundary_ties), expected$boundary_ties)
    scalar <- prime:::indexed_risk_cpp(
      X, y, queries, k, lambda = 1e-6, eps = 0,
      keep_contribution = FALSE, query_block_size = block
    )
    expect_length(scalar$contribution, 0L)
    expect_equal(scalar$risk, actual$risk, tolerance = 1e-12)
    expect_equal(scalar$mean_sse, actual$mean_sse, tolerance = 1e-12)
  }
})

test_that("streaming preserves query ordering and ANN fits across blocks", {
  set.seed(94113)
  X <- matrix(rnorm(111L * 4L), nrow = 111L)
  y <- X[, 1L]^2 + rnorm(nrow(X), sd = 0.2)
  queries <- c(0L, 39L, 110L, 21L, 39L, 7L, 80L)
  permutation <- c(7L, 3L, 5L, 1L, 6L, 2L, 4L)
  for (eps in c(0, 0.3)) {
    actual <- prime:::indexed_risk_cpp(
      X, y, queries, 19L, eps = eps, query_block_size = 1L
    )
    reordered <- prime:::indexed_risk_cpp(
      X, y, queries[permutation], 19L, eps = eps, query_block_size = 4L
    )
    expect_equal(as.numeric(reordered$contribution),
                 as.numeric(actual$contribution)[permutation], tolerance = 1e-12)
    expect_equal(reordered$risk, actual$risk, tolerance = 1e-12)
    expect_equal(actual$contribution[2L], actual$contribution[5L],
                 tolerance = 1e-12)
  }
})

test_that("streamed risk honors boundary ties when responses differ", {
  X <- as.matrix(expand.grid(x1 = -1:1, x2 = -1:1))
  X <- rbind(X, X)
  storage.mode(X) <- "double"
  y <- seq_len(nrow(X))^2 / 100
  queries <- c(0L, 4L, 17L, 9L)
  expected <- indexed_brute_risk(X, y, queries, 7L)
  actual <- prime:::indexed_risk_cpp(X, y, queries, 7L,
                                    query_block_size = 2L)
  expect_true(any(expected$boundary_ties))
  expect_identical(as.logical(actual$boundary_ties), expected$boundary_ties)
  expect_equal(as.numeric(actual$contribution), expected$contribution,
               tolerance = 1e-10)
  expect_equal(actual$risk, expected$risk, tolerance = 1e-10)
})

test_that("forward candidate scores use each projected neighbor geometry", {
  set.seed(94513)
  X <- matrix(rnorm(81L * 5L), nrow = 81L)
  y <- X[, 1L] - X[, 3L]^2 + rnorm(nrow(X), sd = 0.2)
  queries <- c(0L, 41L, 80L, 11L, 41L)
  selected <- c(3L, 2L)
  candidates <- c(5L, 1L, 4L)
  k <- 17L
  expected <- vapply(candidates, function(candidate) {
    columns <- sort(c(selected, candidate))
    risk <- indexed_brute_risk(X[, columns, drop = FALSE], y, queries, k)$risk
    (stats::var(y) - risk) / stats::var(y)
  }, numeric(1))
  for (threads in c(1L, 2L)) {
    actual <- prime:::indexed_forward_values_cpp(
      X, y, queries, selected - 1L, candidates - 1L, k,
      lambda = 1e-6, response_variance = stats::var(y),
      threads = threads, eps = 0, query_block_size = 2L
    )
    expect_equal(as.numeric(actual$values), expected, tolerance = 1e-10)
  }
})

test_that("subset enumeration retains scalar risks and bounded contributions", {
  set.seed(95113)
  X <- matrix(rnorm(65L * 5L), nrow = 65L)
  y <- X[, 1L] - X[, 3L]^2 + rnorm(nrow(X), sd = 0.2)
  data <- prime:::.prime_prepare(X, y, TRUE, NULL, 1e-6,
                                neighbor_eps = 0, query_block_size = 7L)
  k <- 18L
  subsets <- lapply(1:31, function(mask) {
    which(as.logical(intToBits(mask)[seq_len(ncol(X))]))
  })
  for (subset in subsets) {
    expect_true(is.finite(prime:::.prime_risk(data, subset, k)))
  }
  scalar_entries <- as.list(data$risk_cache, all.names = TRUE)
  expect_length(scalar_entries, length(subsets))
  expect_true(all(vapply(scalar_entries, function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value)
  }, logical(1))))
  expect_length(ls(data$value_cache, all.names = TRUE), 0L)

  first <- prime:::.prime_variance_details(data, 1L, k)
  for (subset in subsets[-1L]) {
    details <- prime:::.prime_variance_details(data, subset, k)
    expect_equal(details$risk, prime:::.prime_risk(data, subset, k),
                 tolerance = 1e-12)
    expect_lte(length(ls(data$value_cache, all.names = TRUE)), 2L)
  }
  again <- prime:::.prime_variance_details(data, 1L, k)
  expect_equal(again$contribution, first$contribution, tolerance = 1e-12)
  expect_lte(sum(vapply(as.list(data$value_cache, all.names = TRUE),
                       function(value) length(value$contribution), integer(1))),
             2L * length(data$query_index))
})

test_that("exact default is coherent with explicit bounded query settings", {
  set.seed(5)
  n <- 300L
  X <- matrix(runif(n * 8L, -1, 1), nrow = n)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  y <- 2 * X[, 1L] - X[, 2L] + 1.5 * X[, 4L]^2 + rnorm(n, sd = 0.1)
  default <- prime(X, y)
  explicit <- prime(X, y, neighbor_eps = 0, query_block_size = 7L)
  expect_identical(default$selected, explicit$selected)
  expect_equal(default$effects, explicit$effects, tolerance = 1e-10)
  expect_equal(default$subset_values, explicit$subset_values, tolerance = 1e-10)
  expect_equal(sum(default$effects), 1, tolerance = 1e-10)
  expect_true(all(c(1L, 2L, 4L) %in% default$selected))
  expect_identical(default$settings$neighbor_eps, 0)
  expect_identical(explicit$settings$query_block_size, 7L)
})

test_that("public ANN and block settings reject invalid values clearly", {
  set.seed(96113)
  X <- matrix(rnorm(30L * 2L), nrow = 30L)
  y <- X[, 1L] + rnorm(30L)
  invalid_eps <- list(-0.1, Inf, NA_real_, c(0, 0.1), numeric(), "0.1", TRUE)
  invalid_blocks <- list(0, -1, 1.5, Inf, NA_real_, c(1L, 2L),
                         numeric(), "2", TRUE, 2^31)
  for (function_name in c("prime", "prime_screen")) {
    fit_function <- get(function_name, asNamespace("prime"))
    for (value in invalid_eps) {
      expect_error(fit_function(X, y, neighbor_eps = value), "neighbor_eps")
    }
    for (value in invalid_blocks) {
      expect_error(fit_function(X, y, query_block_size = value),
                   "query_block_size")
    }
  }
})
