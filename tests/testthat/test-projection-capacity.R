test_that("projection capacity is checked arithmetically without allocating edges", {
  count <- prime:::.prime_projection_edge_count
  expect_identical(count(27L), 27 * 2^26)
  expect_lte(count(27L), .Machine$integer.max)
  for (dimension in 28:30) {
    expect_error(count(dimension), "projection edge count exceeds 32-bit integer-index capacity", fixed = TRUE)
  }
  expect_identical(prime:::.prime_max_mask_dimension, 30L)
})

test_that("public allocation rejects excessive edges before input materialization", {
  allocate <- prime_allocate
  scope <- new.env(parent = environment(allocate))
  scope$length <- function(x) {
    if (inherits(x, "virtual_subset_values")) 2^28 else base::length(x)
  }
  scope$as.numeric <- function(...) stop("subset input was materialized")
  environment(allocate) <- scope
  virtual <- structure(numeric(), class = "virtual_subset_values")
  expect_error(allocate(virtual), "projection edge count exceeds 32-bit integer-index capacity", fixed = TRUE)
})

test_that("prime rejects excessive edges before subset evaluation", {
  fit <- prime
  scope <- new.env(parent = environment(fit))
  scope$prime_screen <- function(...) {
    list(selected = 1:28, feature_names = paste0("X", 1:28))
  }
  scope$.prime_estimate_subset_values <- function(...) {
    stop("subset enumeration was reached")
  }
  environment(fit) <- scope
  expect_error(fit(NULL, NULL), "projection edge count exceeds 32-bit integer-index capacity", fixed = TRUE)
})

test_that("internal constructor and projection reject before evaluating data", {
  expect_error(
    prime:::.prime_estimate_subset_values(stop("risk data was evaluated"), 1:28),
    "projection edge count exceeds 32-bit integer-index capacity", fixed = TRUE
  )
  expect_error(
    prime:::.prime_project_monotone(stop("game was evaluated"), 28L),
    "projection edge count exceeds 32-bit integer-index capacity", fixed = TRUE
  )
  edges <- prime:::.prime_hasse_edges
  scope <- new.env(parent = environment(edges))
  scope$integer <- function(...) stop("edge storage was allocated")
  environment(edges) <- scope
  expect_error(edges(28L), "projection edge count exceeds 32-bit integer-index capacity", fixed = TRUE)
})
