.prime_empty_forward_trace <- function() {
  data.frame(
    step = integer(),
    added = integer(),
    added_name = character(),
    selected_before = character(),
    dimension = integer(),
    K = integer(),
    query_count = integer(),
    current_value = double(),
    candidate_value = double(),
    gain = double(),
    raw_gain = double(),
    query_scale = double(),
    threshold = double(),
    calibrated_gain = double(),
    accepted = logical(),
    stringsAsFactors = FALSE
  )
}

.prime_empty_backward_trace <- function() {
  data.frame(
    step = integer(),
    tested = integer(),
    tested_name = character(),
    selected_before = character(),
    dimension = integer(),
    K = integer(),
    query_count = integer(),
    full_risk = double(),
    raw_loss = double(),
    clipped_loss = double(),
    query_scale = double(),
    threshold = double(),
    calibrated_loss = double(),
    removed = logical(),
    stringsAsFactors = FALSE
  )
}

.prime_ridge_initial <- function(data) {
  ridge_cap <- .prime_ridge_cap(data$n, data$p)
  # Only the ridge/BIC calculation uses response units standardized to variance 1.
  # This keeps its numerical floor and coefficient ranking invariant to units.
  scaled_y <- (data$y - mean(data$y)) / sqrt(data$response_variance)
  frame <- data.frame(y = scaled_y, data$X, check.names = FALSE)
  names(frame) <- c("y", paste0("X", seq_len(data$p)))
  fit <- if (data$p == 1L) {
    centered_x <- data$X[, 1L] - mean(data$X[, 1L])
    scaled_x <- centered_x / sqrt(mean(centered_x^2))
    centered_y <- scaled_y - mean(scaled_y)
    # MASS uses centered, root-mean-square-scaled predictor coordinates.
    # Its intercept path drops the matrix dimension for a single predictor.
    list(coef = sum(scaled_x * centered_y) /
      (sum(scaled_x^2) + data$ridge_lambda))
  } else {
    MASS::lm.ridge(
      y ~ ., data = frame, lambda = data$ridge_lambda
    )
  }
  beta_squared <- as.numeric(fit$coef)^2
  total <- sum(beta_squared)
  weights <- if (is.finite(total) && total > 0) {
    beta_squared / total
  } else {
    rep(0, data$p)
  }
  ordered <- order(-weights, seq_along(weights))
  path <- lapply(0:ridge_cap, function(size) {
    if (!size) {
      residual <- scaled_y - mean(scaled_y)
      rank <- 0L
    } else {
      fit <- stats::lm.fit(
        x = cbind(1, data$X[, ordered[seq_len(size)], drop = FALSE]),
        y = scaled_y
      )
      residual <- fit$residuals
      rank <- max(0L, fit$rank - 1L)
    }
    scaled_rss <- sum(residual^2)
    data.frame(
      size = size,
      rank = rank,
      rss = scaled_rss * data$response_variance,
      bic = data$n * log(max(scaled_rss / data$n, .Machine$double.eps)) +
        rank * log(data$n)
    )
  })
  path <- do.call(rbind, path)
  chosen_size <- path$size[which.min(path$bic)]
  selected <- if (chosen_size) {
    sort(ordered[seq_len(chosen_size)])
  } else {
    integer()
  }
  list(
    selected = as.integer(selected),
    weights = weights,
    cap = ridge_cap,
    chosen_size = chosen_size,
    path = path
  )
}

.prime_forward_core <- function(data) {
  ridge <- .prime_ridge_initial(data)
  selected <- sort(unique(as.integer(ridge$selected)))
  trace <- list()
  dimension_limit <- .prime_dimension_limit(data$n, data$p)
  stop_reason <- "all coordinates selected"

  repeat {
    if (length(selected) >= data$p) break
    if (length(selected) >= dimension_limit) {
      stop_reason <- "local-fit dimension limit reached"
      break
    }
    candidates <- setdiff(seq_len(data$p), selected)
    dimension <- length(selected) + 1L
    k <- .prime_neighborhood_size(data$n, dimension)
    if (!is.finite(k)) {
      stop_reason <- "local-fit dimension limit reached"
      break
    }

    current_value <- if (length(selected)) {
      (data$response_variance - .prime_variance_details(data, selected, k)$risk) /
        data$response_variance
    } else {
      0
    }
    candidate_values <- .prime_forward_values(
      data, selected, candidates, k
    )
    best_position <- which.max(candidate_values)
    best_value <- candidate_values[best_position]
    gain <- best_value - current_value
    if (!length(best_position) || !is.finite(best_value)) {
      stop("No finite forward candidate value was available.",
           call. = FALSE)
    }

    added <- candidates[best_position]
    diagnostic <- .prime_forward_diagnostic(
      data, selected, added, k
    )
    accepted <- is.finite(diagnostic$normalized_calibrated_gain) &&
      diagnostic$normalized_calibrated_gain > 0
    trace[[length(trace) + 1L]] <- data.frame(
      step = length(trace) + 1L,
      added = added,
      added_name = data$feature_names[added],
      selected_before = paste(selected, collapse = ","),
      dimension = dimension,
      K = k,
      query_count = length(data$query_index),
      current_value = current_value,
      candidate_value = best_value,
      gain = gain,
      raw_gain = diagnostic$raw_gain,
      query_scale = diagnostic$query_scale,
      threshold = diagnostic$threshold,
      calibrated_gain = diagnostic$calibrated_gain,
      accepted = accepted,
      stringsAsFactors = FALSE
    )
    if (!accepted) {
      stop_reason <- "raw-gain leader's calibrated gain was nonpositive"
      break
    }
    selected <- sort(c(selected, added))
  }

  candidates <- as.integer(selected)
  local_path <- setdiff(candidates, ridge$selected)
  list(
    candidates = as.integer(candidates),
    candidate_names = data$feature_names[candidates],
    ridge_initial = ridge$selected,
    ridge_initial_names = data$feature_names[ridge$selected],
    ridge_weights = stats::setNames(ridge$weights, data$feature_names),
    ridge_cap = ridge$cap,
    ridge_chosen_size = ridge$chosen_size,
    ridge_path = ridge$path,
    ridge_lambda = data$ridge_lambda,
    local_path = local_path,
    local_path_names = data$feature_names[local_path],
    dimension_limit = dimension_limit,
    additions = length(local_path),
    trace = if (length(trace)) {
      do.call(rbind, trace)
    } else {
      .prime_empty_forward_trace()
    },
    stop_reason = stop_reason
  )
}

.prime_backward_core <- function(data, candidates) {
  selected <- sort(unique(as.integer(candidates)))
  trace <- list()
  if (!length(selected)) {
    return(list(
      selected = integer(),
      selected_names = character(),
      trace = .prime_empty_backward_trace()
    ))
  }

  repeat {
    d <- length(selected)
    k <- .prime_neighborhood_size(data$n, d)
    if (!is.finite(k)) {
      stop("Backward pruning is infeasible at dimension ", d, ".",
           call. = FALSE)
    }
    full <- .prime_variance_details(data, selected, k)
    diagnostics <- lapply(selected, function(variable) {
      diagnostic <- .prime_deletion_diagnostic(
        data, selected, variable, k, full_details = full
      )
      diagnostic$paired_loss <- NULL
      diagnostic
    })
    calibrated <- vapply(
      diagnostics, function(value) value$normalized_calibrated_loss, numeric(1)
    )
    tested_position <- which.min(calibrated)
    tested <- selected[tested_position]
    tested_diagnostic <- diagnostics[[tested_position]]
    remove <- is.finite(calibrated[tested_position]) &&
      calibrated[tested_position] == 0
    trace[[length(trace) + 1L]] <- data.frame(
      step = length(trace) + 1L,
      tested = tested,
      tested_name = data$feature_names[tested],
      selected_before = paste(selected, collapse = ","),
      dimension = d,
      K = k,
      query_count = length(data$query_index),
      full_risk = full$risk,
      raw_loss = tested_diagnostic$raw_loss,
      clipped_loss = tested_diagnostic$clipped_loss,
      query_scale = tested_diagnostic$query_scale,
      threshold = tested_diagnostic$threshold,
      calibrated_loss = tested_diagnostic$calibrated_loss,
      removed = remove,
      stringsAsFactors = FALSE
    )
    if (!remove) break
    selected <- setdiff(selected, tested)
    if (!length(selected)) break
  }

  list(
    selected = as.integer(selected),
    selected_names = data$feature_names[selected],
    trace = if (length(trace)) {
      do.call(rbind, trace)
    } else {
      .prime_empty_backward_trace()
    }
  )
}

#' Screen coordinates by local-linear explained variance
#'
#' @param X Numeric predictor matrix or data frame.
#' @param y Numeric response.
#' @param rescale Standardize predictors before distance calculations.
#' @param twin_r Inverse Data-Twinning split ratio; \code{NULL} uses all
#'   observations as query locations.
#' @param lambda Numerical slope-only nugget used in local fits.
#' @param neighbor_eps Nonnegative relative Euclidean neighbor-distance tolerance;
#'   zero uses exact indexed search and positive values enable approximation.
#' @param query_block_size Positive query batch size; the native risk engine
#'   streams neighborhoods and checks interrupts at these intervals.
#' @return An object containing forward candidates, selected coordinates,
#'   diagnostic traces, and fixed settings.
#' @export
prime_screen <- function(
  X,
  y,
  rescale = TRUE,
  twin_r = NULL,
  lambda = 1e-6,
  neighbor_eps = 0,
  query_block_size = 256L
) {
  call <- match.call()
  data <- .prime_prepare(X, y, rescale, twin_r, lambda,
                         neighbor_eps, query_block_size)
  forward <- .prime_forward_core(data)
  backward <- .prime_backward_core(data, forward$candidates)
  structure(
    list(
      call = call,
      candidates = forward$candidates,
      candidate_names = forward$candidate_names,
      selected = backward$selected,
      selected_names = backward$selected_names,
      forward = forward,
      backward = backward,
      settings = .prime_settings(rescale, twin_r, lambda,
                                 neighbor_eps, query_block_size),
      feature_names = data$feature_names,
      data = data
    ),
    class = "prime_screen"
  )
}

#' @export
print.prime_screen <- function(x, ...) {
  retained <- paste(x$selected_names, collapse = ", ")
  if (!nzchar(retained)) retained <- "<none>"
  cat("PRIME screen\n")
  cat("  Candidates:", length(x$candidate_names), "\n")
  cat("  Retained:", length(x$selected_names), "(", retained, ")\n")
  cat("  Use summary() for forward and backward diagnostics.\n")
  invisible(x)
}

#' Summarize a PRIME screening result
#'
#' @param object A \code{prime_screen} object.
#' @param ... Unused.
#' @return A structured summary of the forward and backward selection stages.
#' @method summary prime_screen
#' @export
summary.prime_screen <- function(object, ...) {
  feature_names <- object$feature_names
  to_names <- function(index) {
    if (length(index)) feature_names[index] else character()
  }
  backward_trace <- object$backward$trace
  structure(
    list(
      call = object$call,
      n = object$data$n,
      p = object$data$p,
      ridge_seed = to_names(object$forward$ridge_initial),
      local_path = to_names(object$forward$local_path),
      candidates = object$candidate_names,
      retained = object$selected_names,
      forward_stop = object$forward$stop_reason,
      backward_evaluations = if (is.data.frame(backward_trace)) {
        nrow(backward_trace)
      } else {
        length(backward_trace)
      },
      backward_deletions = if (
        is.data.frame(backward_trace) && "removed" %in% names(backward_trace)
      ) {
        as.integer(sum(backward_trace$removed))
      } else {
        0L
      },
      settings = object$settings
    ),
    class = "summary_prime_screen"
  )
}

#' @export
print.summary_prime_screen <- function(x, ...) {
  format_names <- function(values) {
    text <- paste(values, collapse = ", ")
    if (nzchar(text)) text else "<none>"
  }
  cat("Call:\n")
  print(x$call)
  cat("\nData:\n")
  cat(" ", x$n, " observations, ", x$p, " predictors\n", sep = "")
  cat("\nForward selection:\n")
  cat(" Ridge seed:", format_names(x$ridge_seed), "\n")
  cat(" Local path:", format_names(x$local_path), "\n")
  cat(" Stop:", x$forward_stop, "\n")
  cat("\nBackward reduction:\n")
  cat(" Candidates:", format_names(x$candidates), "\n")
  cat(" Retained:", format_names(x$retained), "\n")
  cat(" Evaluations:", x$backward_evaluations, "\n")
  cat(" Deletions:", x$backward_deletions, "\n")
  invisible(x)
}
