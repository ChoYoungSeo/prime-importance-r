.prime_validate_twin_r <- function(twin_r, n) {
  if (is.null(twin_r)) return(NULL)
  value <- suppressWarnings(as.numeric(twin_r))
  if (length(value) != 1L || !is.finite(value) ||
      value != floor(value) || value < 2L || value > floor(n / 2L)) {
    stop(
      "twin_r must be NULL or an integer between 2 and floor(n/2).",
      call. = FALSE
    )
  }
  as.integer(value)
}

.prime_prepare <- function(X, y, rescale, twin_r, lambda,
                           neighbor_eps = 0, query_block_size = 256L) {
  if (!is.numeric(neighbor_eps) || length(neighbor_eps) != 1L ||
      !is.finite(neighbor_eps) || neighbor_eps < 0 || neighbor_eps > 1e18) {
    stop("neighbor_eps must be a finite nonnegative number at most 1e18.",
         call. = FALSE)
  }
  if (!is.numeric(query_block_size) || length(query_block_size) != 1L ||
      !is.finite(query_block_size) || query_block_size < 1 ||
      query_block_size != floor(query_block_size) ||
      query_block_size > .Machine$integer.max) {
    stop("query_block_size must be a positive integer within the R integer range.",
         call. = FALSE)
  }
  query_block_size <- as.integer(query_block_size)

  if (is.data.frame(X) && !all(vapply(X, is.numeric, logical(1)))) {
    stop("All columns of X must be numeric.", call. = FALSE)
  }
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  y <- as.numeric(y)

  if (!is.logical(rescale) || length(rescale) != 1L || is.na(rescale)) {
    stop("rescale must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(dim(X)) != 2L || nrow(X) < 10L || ncol(X) < 1L) {
    stop("X must have at least 10 rows and one column.", call. = FALSE)
  }
  if (nrow(X) != length(y)) {
    stop("nrow(X) must equal length(y).", call. = FALSE)
  }
  if (any(!is.finite(X)) || any(!is.finite(y))) {
    stop("X and y must contain only finite values.", call. = FALSE)
  }
  predictor_scale <- apply(X, 2L, stats::sd)
  if (any(!is.finite(predictor_scale) | predictor_scale == 0)) {
    stop("X must not contain constant columns.", call. = FALSE)
  }
  response_variance <- stats::var(y)
  if (!is.finite(response_variance) || response_variance <= 0) {
    stop("y must be nonconstant.", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda < 0) {
    stop("lambda must be a finite nonnegative scalar.", call. = FALSE)
  }

  feature_names <- colnames(X)
  if (is.null(feature_names)) {
    feature_names <- paste0("X", seq_len(ncol(X)))
  }
  if (anyDuplicated(feature_names)) {
    stop("Column names of X must be unique.", call. = FALSE)
  }

  center <- rep(0, ncol(X))
  column_scale <- rep(1, ncol(X))
  if (rescale) {
    center <- colMeans(X)
    column_scale <- predictor_scale
    X <- sweep(sweep(X, 2L, center, "-"), 2L, column_scale, "/")
  }

  twin_r <- .prime_validate_twin_r(twin_r, nrow(X))
  query_index <- seq_len(nrow(X))
  if (!is.null(twin_r)) {
    if (!requireNamespace("twinning", quietly = TRUE)) {
      stop(
        "Package 'twinning' is required when twin_r is specified.",
        call. = FALSE
      )
    }
    query_design <- cbind(X, y = y)
    query_scale <- apply(query_design, 2L, stats::sd)
    query_design <- query_design[
      ,
      is.finite(query_scale) & query_scale > 0,
      drop = FALSE
    ]
    query_design <- scale(query_design)
    query_index <- as.integer(twinning::twin(
      query_design,
      r = twin_r,
      u1 = which.max(rowSums(query_design^2)),
      format_data = FALSE
    ))
  }

  list(
    X = X,
    y = y,
    n = nrow(X),
    p = ncol(X),
    response_variance = response_variance,
    feature_names = feature_names,
    center = center,
    scale = column_scale,
    query_index = query_index,
    twin_r = twin_r,
    value_cache = new.env(parent = emptyenv()),
    risk_cache = new.env(parent = emptyenv()),
    cache_state = list2env(list(order = character()), parent = emptyenv()),
    neighbor_eps = neighbor_eps,
    query_block_size = query_block_size,
    rescale = rescale,
    ridge_lambda = .prime_constants$ridge_lambda,
    lambda = as.numeric(lambda)
  )
}
