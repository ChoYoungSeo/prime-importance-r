all_permutations <- function(values) {
  if (length(values) == 1L) return(matrix(values, nrow = 1L))
  do.call(rbind, lapply(seq_along(values), function(i) {
    cbind(values[i], all_permutations(values[-i]))
  }))
}

permutation_pv <- function(game, s) {
  permutations <- all_permutations(seq_len(s))
  allocations <- matrix(0, nrow(permutations), s)
  log_weight <- numeric(nrow(permutations))
  for (row in seq_len(nrow(permutations))) {
    mask <- 0L
    previous <- game[1L]
    for (step in seq_len(s)) {
      variable <- permutations[row, step]
      mask <- bitwOr(mask, bitwShiftL(1L, variable - 1L))
      current <- game[mask + 1L]
      allocations[row, variable] <- current - previous
      log_weight[row] <- log_weight[row] - log(current)
      previous <- current
    }
  }
  weight <- exp(log_weight - max(log_weight))
  colSums(allocations * weight) / sum(weight)
}

test_that("ratio-potential PME equals complete random-order enumeration", {
  s <- 5L
  weights <- c(0.4, 0.7, 1.1, 1.6, 2.0)
  masks <- 0:(2^s - 1L)
  game <- vapply(masks, function(mask) {
    active <- vapply(
      seq_len(s),
      function(j) bitwAnd(mask, bitwShiftL(1L, j - 1L)) != 0L,
      logical(1)
    )
    if (!any(active)) return(0)
    sum(weights[active])^1.3
  }, numeric(1))
  expected <- permutation_pv(game, s)
  actual <- unclass(prime_allocate(game))
  attributes(actual) <- NULL

  expect_equal(actual, expected, tolerance = 1e-11)
  expect_equal(sum(actual), game[length(game)], tolerance = 1e-11)
})

test_that("additive games return their additive weights", {
  weights <- c(0.1, 0.2, 0.3, 0.4)
  s <- length(weights)
  game <- vapply(0:(2^s - 1L), function(mask) {
    sum(weights[vapply(
      seq_len(s),
      function(j) bitwAnd(mask, bitwShiftL(1L, j - 1L)) != 0L,
      logical(1)
    )])
  }, numeric(1))
  actual <- prime_allocate(game)
  expect_equal(as.numeric(actual), weights, tolerance = 1e-12)
})

test_that("PME excludes a correlated exogenous proxy", {
  rho <- 0.9
  game <- c(
    0,
    (1 - rho^2) / 2,
    1 / 2,
    1 - rho^2 / 2,
    0,
    1 / 2,
    1 / 2,
    1
  )
  actual <- prime_allocate(game)
  expect_equal(as.numeric(actual), c(0.5, 0.5, 0), tolerance = 1e-12)
})

test_that("PME on a sufficient superset equals PME on the active game", {
  active_weights <- c(0.2, 0.5, 0.8, 1.1)
  active_size <- length(active_weights)
  proxy_size <- 3L
  active_game <- vapply(0:(2^active_size - 1L), function(mask) {
    present <- vapply(seq_len(active_size), function(position) {
      bitwAnd(mask, bitwShiftL(1L, position - 1L)) != 0L
    }, logical(1))
    if (!any(present)) return(0)
    sum(active_weights[present])^1.2
  }, numeric(1))
  active_game <- active_game / active_game[length(active_game)]
  full_size <- active_size + proxy_size
  full_game <- vapply(0:(2^full_size - 1L), function(mask) {
    active_mask <- bitwAnd(mask, 2^active_size - 1L)
    active_game[active_mask + 1L]
  }, numeric(1))

  reduced <- prime_allocate(active_game)
  full <- prime_allocate(full_game)
  expect_equal(
    as.numeric(full[seq_len(active_size)]),
    as.numeric(reduced),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(full[active_size + seq_len(proxy_size)]),
    numeric(proxy_size),
    tolerance = 1e-12
  )
})

test_that("monotone projection respects every Boolean-lattice edge", {
  raw <- c(0, 1.2, -0.1, 0.9, 0.4, 0.8, 1.1, 1)
  effects <- prime_allocate(raw)
  projected <- attr(effects, "projected_values")
  edges <- prime:::.prime_hasse_edges(3L)
  expect_gte(
    min(projected[edges$upper + 1L] - projected[edges$lower + 1L]),
    -1e-9
  )
  expect_equal(projected[1L], 0)
  expect_equal(projected[length(projected)], 1)
  expect_equal(sum(effects), 1, tolerance = 1e-9)
})
