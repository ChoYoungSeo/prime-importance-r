.prime_max_mask_dimension <- 30L

.prime_constants <- list(
  ridge_cap_multiplier = 1.5,
  neighborhood_a = 1.5,
  neighborhood_b = 3,
  cap_fraction = 0.5,
  ridge_lambda = 1e-6
)

.prime_settings <- function(rescale, twin_r, lambda, neighbor_eps = 0,
                            query_block_size = 256L) {
  query_rule <- if (is.null(twin_r)) {
    "all observed query locations"
  } else {
    paste0(
      "Data Twinning on standardized (X,y), inverse split ratio r=",
      twin_r
    )
  }
  list(
    ridge_rule = paste0(
      "BIC-selected size along the standardized ridge-weight ranking, ",
      "searched through ceiling(1.5*log(max(2,p))) coordinates"
    ),
    ridge_lambda = .prime_constants$ridge_lambda,
    forward_rule = paste0(
      "start from the ridge seed and extend the local-linear path while the ",
      "leading query-calibrated gain is positive, stopped at the local-fit ",
      "dimension limit"
    ),
    backward_rule = paste0(
      "sequentially delete zero calibrated losses using paired querywise ",
      "empirical scales and the sqrt(2*log(max(2,dimension))) universal threshold"
    ),
    neighborhood_rule = "ceiling(max(1.5*sqrt(n), 3*d))",
    neighborhood_cap = "floor(n/2)",
    neighbor_search = if (neighbor_eps == 0) "exact projected nearest neighbors"
      else "approximate projected nearest neighbors",
    neighbor_eps = neighbor_eps,
    query_block_size = as.integer(query_block_size),
    neighbor_ties = if (neighbor_eps == 0) "distance then original row index"
      else "returned neighbors ordered by distance then original row index",
    neighbor_backend = "nanoflann 1.9.0 indexed streaming",
    query_rule = query_rule,
    twin_r = twin_r,
    local_lambda = lambda,
    lambda = lambda,
    residual_df = "K-d-1",
    rescale = rescale
  )
}

.prime_dimension_limit <- function(n, p) {
  as.integer(min(
    p,
    floor(n * .prime_constants$cap_fraction /
            .prime_constants$neighborhood_b)
  ))
}

.prime_ridge_cap <- function(n, p) {
  target <- max(
    1L,
    ceiling(.prime_constants$ridge_cap_multiplier * log(max(2L, p)))
  )
  as.integer(min(p, target, .prime_dimension_limit(n, p)))
}

.prime_neighborhood_size <- function(n, d) {
  k <- ceiling(max(
    .prime_constants$neighborhood_a * sqrt(n),
    .prime_constants$neighborhood_b * d
  ))
  cap <- min(n - 1L, floor(.prime_constants$cap_fraction * n))
  if (!is.finite(k) || k <= d + 1L || k > cap) {
    return(NA_integer_)
  }
  as.integer(k)
}

.prime_worker_count <- function(candidate_count) {
  configured <- suppressWarnings(as.integer(
    Sys.getenv("PRIME_NUM_THREADS", "")
  ))
  if (length(configured) == 1L && is.finite(configured) &&
      configured >= 1L) {
    return(as.integer(max(1L, min(candidate_count, configured, 8L))))
  }
  cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!length(cores) || !is.finite(cores) || cores < 1L) cores <- 1L
  as.integer(max(1L, min(candidate_count, cores, 8L)))
}
