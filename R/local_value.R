.prime_query_neighbors <- function(X, query_index, k) {
  algorithm <- if (nrow(X) <= 2000L && ncol(X) >= 5L) {
    "brute"
  } else {
    "kd_tree"
  }
  all_queries <- length(query_index) == nrow(X) &&
    identical(query_index, seq_len(nrow(X)))
  if (all_queries) {
    return(FNN::get.knn(
      X,
      k = as.integer(k),
      algorithm = algorithm
    )$nn.index)
  }

  candidates <- FNN::get.knnx(
    X,
    X[query_index, , drop = FALSE],
    k = as.integer(k) + 1L,
    algorithm = algorithm
  )$nn.index
  query_matrix <- matrix(
    query_index,
    nrow = length(query_index),
    ncol = ncol(candidates)
  )
  self_match <- candidates == query_matrix
  self_position <- max.col(self_match, ties.method = "first")
  self_position[rowSums(self_match) == 0L] <- ncol(candidates) + 1L

  target_column <- matrix(
    rep(seq_len(k), each = length(query_index)),
    nrow = length(query_index),
    ncol = k
  )
  source_column <- target_column + (target_column >= self_position)
  row_index <- matrix(
    rep(seq_along(query_index), times = k),
    nrow = length(query_index),
    ncol = k
  )
  neighbors <- matrix(
    candidates[cbind(as.vector(row_index), as.vector(source_column))],
    nrow = length(query_index),
    ncol = k
  )
  if (any(neighbors == matrix(
    query_index,
    nrow = length(query_index),
    ncol = k
  ))) {
    stop("A query observation remained in its own neighborhood.",
         call. = FALSE)
  }
  neighbors
}

.prime_compute_indexed <- function(data, subset, k, keep_contribution) {
  subset <- sort(unique(as.integer(subset)))
  if (!length(subset)) {
    return(list(
      risk = data$response_variance, mean_sse = NA_real_, residual_df = NA_real_,
      contribution = if (keep_contribution) rep(data$response_variance,
        length(data$query_index)) else numeric(),
      query_count = length(data$query_index)
    ))
  }
  details <- indexed_risk_cpp(
    data$X[, subset, drop = FALSE], data$y, data$query_index - 1L,
    as.integer(k), data$lambda, data$neighbor_eps,
    keep_contribution, data$query_block_size
  )
  details$contribution <- as.numeric(details$contribution)
  details$query_count <- length(data$query_index)
  details
}

.prime_variance_details <- function(data, subset, k) {
  subset <- sort(unique(as.integer(subset)))
  key <- paste(k, paste(subset, collapse = ","), sep = ":")
  if (exists(key, envir = data$value_cache, inherits = FALSE)) {
    data$cache_state$order <- c(setdiff(data$cache_state$order, key), key)
    return(get(key, envir = data$value_cache, inherits = FALSE))
  }
  details <- .prime_compute_indexed(data, subset, k, TRUE)
  assign(key, details$risk, envir = data$risk_cache)
  assign(key, details, envir = data$value_cache)
  data$cache_state$order <- c(setdiff(data$cache_state$order, key), key)
  if (length(data$cache_state$order) > 2L) {
    expired <- utils::head(data$cache_state$order, -2L)
    rm(list = expired, envir = data$value_cache)
    data$cache_state$order <- utils::tail(data$cache_state$order, 2L)
  }
  details
}

.prime_risk <- function(data, subset, k) {
  subset <- sort(unique(as.integer(subset)))
  if (!length(subset)) return(data$response_variance)
  key <- paste(k, paste(subset, collapse = ","), sep = ":")
  if (exists(key, envir = data$risk_cache, inherits = FALSE)) {
    return(get(key, envir = data$risk_cache, inherits = FALSE))
  }
  risk <- .prime_compute_indexed(data, subset, k, FALSE)$risk
  assign(key, risk, envir = data$risk_cache)
  risk
}

.prime_value <- function(data, subset, k) {
  if (!length(subset)) return(0)
  risk <- .prime_risk(data, subset, k)
  max(0, data$response_variance - risk) / data$response_variance
}

.prime_forward_score <- function(data, subset, k) {
  if (!length(subset)) return(0)
  risk <- .prime_risk(data, subset, k)
  (data$response_variance - risk) / data$response_variance
}

.prime_deletion_diagnostic <- function(
  data, selected, variable, k, full_details = NULL
) {
  selected <- sort(unique(as.integer(selected)))
  if (!variable %in% selected) {
    stop("The tested coordinate is not in the selected set.", call. = FALSE)
  }
  if (is.null(full_details)) {
    full_details <- .prime_variance_details(data, selected, k)
  }
  reduced <- setdiff(selected, variable)
  reduced_details <- .prime_variance_details(data, reduced, k)
  paired_loss <- reduced_details$contribution - full_details$contribution
  scaled_loss <- paired_loss / data$response_variance
  raw_scaled <- mean(scaled_loss)
  scale_scaled <- if (length(scaled_loss) > 1L) {
    stats::sd(scaled_loss) / sqrt(length(scaled_loss))
  } else {
    0
  }
  threshold_scaled <- scale_scaled * sqrt(
    2 * log(max(2L, length(selected)))
  )
  raw_loss <- raw_scaled * data$response_variance
  query_scale <- scale_scaled * data$response_variance
  threshold <- threshold_scaled * data$response_variance
  list(
    raw_loss = raw_loss,
    clipped_loss = max(raw_loss, 0),
    query_scale = query_scale,
    threshold = threshold,
    calibrated_loss = max(raw_scaled - threshold_scaled, 0) *
      data$response_variance,
    normalized_calibrated_loss = max(raw_scaled - threshold_scaled, 0),
    full_risk = full_details$risk,
    reduced_risk = reduced_details$risk,
    paired_loss = paired_loss
  )
}

.prime_forward_diagnostic <- function(
  data, selected, variable, k
) {
  selected <- sort(unique(as.integer(selected)))
  current <- .prime_variance_details(data, selected, k)
  augmented <- .prime_variance_details(
    data, sort(c(selected, variable)), k
  )
  paired_gain <- current$contribution - augmented$contribution
  scaled_gain <- paired_gain / data$response_variance
  raw_scaled <- mean(scaled_gain)
  scale_scaled <- if (length(scaled_gain) > 1L) {
    stats::sd(scaled_gain) / sqrt(length(scaled_gain))
  } else {
    0
  }
  augmented_dimension <- length(selected) + 1L
  threshold_scaled <- scale_scaled * sqrt(
    2 * log(max(2L, augmented_dimension))
  )
  raw_gain <- raw_scaled * data$response_variance
  query_scale <- scale_scaled * data$response_variance
  threshold <- threshold_scaled * data$response_variance
  list(
    raw_gain = raw_gain,
    query_scale = query_scale,
    threshold = threshold,
    calibrated_gain = max(raw_scaled - threshold_scaled, 0) *
      data$response_variance,
    normalized_calibrated_gain = max(raw_scaled - threshold_scaled, 0),
    paired_gain = paired_gain
  )
}

.prime_forward_values <- function(data, selected, candidates, k) {
  batch <- indexed_forward_values_cpp(
    data$X, data$y, data$query_index - 1L,
    as.integer(selected) - 1L, as.integer(candidates) - 1L,
    as.integer(k), data$lambda, data$response_variance,
    .prime_worker_count(length(candidates)), data$neighbor_eps,
    data$query_block_size
  )
  as.numeric(batch$values)
}

.prime_variance_details_r <- function(data, subset, k) {
  subset <- sort(unique(as.integer(subset)))
  if (!length(subset)) {
    return(list(
      risk = data$response_variance,
      mean_sse = NA_real_,
      residual_df = NA_real_,
      contribution = rep(
        data$response_variance, length(data$query_index)
      )
    ))
  }
  Xs <- data$X[, subset, drop = FALSE]
  neighbors <- .prime_query_neighbors(Xs, data$query_index, k)
  local_sse <- vapply(seq_along(data$query_index), function(i) {
    index <- neighbors[i, ]
    centered <- sweep(
      Xs[index, , drop = FALSE],
      2L,
      Xs[data$query_index[i], ],
      "-"
    )
    design <- cbind(1, centered)
    cross <- crossprod(design)
    diag(cross)[-1L] <- diag(cross)[-1L] + data$lambda
    coefficient <- tryCatch(
      solve(cross, crossprod(design, data$y[index])),
      error = function(e) {
        qr.solve(cross, crossprod(design, data$y[index]))
      }
    )
    sum((data$y[index] - design %*% coefficient)^2)
  }, numeric(1))
  residual_df <- k - length(subset) - 1
  mean_sse <- mean(local_sse)
  list(
    risk = mean_sse / residual_df,
    mean_sse = mean_sse,
    residual_df = residual_df,
    contribution = local_sse / residual_df
  )
}
