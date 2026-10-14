#' Proportional inclusive marginal effects
#'
#' Exact subset enumeration scales exponentially with the selected dimension.
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
#' @return A fitted screen, reduced PRIME allocation, and subset-value
#'   diagnostics. If screening retains no coordinate, \code{effects} is a
#'   zero vector and \code{subset_values} is \code{NULL}, indicating that a
#'   normalized allocation is unavailable.
#' @export
prime <- function(
  X,
  y,
  rescale = TRUE,
  twin_r = NULL,
  lambda = 1e-6,
  neighbor_eps = 0,
  query_block_size = 256L
) {
  call <- match.call()
  screen <- prime_screen(
    X, y,
    rescale = rescale,
    twin_r = twin_r,
    lambda = lambda,
    neighbor_eps = neighbor_eps,
    query_block_size = query_block_size
  )
  screen$call <- call
  selected <- screen$selected
  full_effects <- stats::setNames(
    numeric(length(screen$feature_names)),
    screen$feature_names
  )
  if (!length(selected)) {
    return(structure(
      list(
        call = call,
        effects = full_effects,
        selected = integer(),
        selected_names = character(),
        screen = screen,
        subset_values = NULL,
        signal = 0,
        noise_floor = screen$data$response_variance,
        risk_full = screen$data$response_variance,
        projection = NULL,
        settings = screen$settings
      ),
      class = "prime_fit"
    ))
  }

  .prime_projection_edge_count(length(selected))
  # Allocation consumes scalar risks, not the selection contribution cache.
  rm(list = ls(screen$data$value_cache, all.names = TRUE),
     envir = screen$data$value_cache)
  screen$data$cache_state$order <- character()
  estimated_values <- .prime_estimate_subset_values(screen$data, selected)
  reduced_effects <- prime_allocate(estimated_values$raw)
  names(reduced_effects) <- screen$feature_names[selected]
  full_effects[selected] <- reduced_effects
  projected <- attr(reduced_effects, "projected_values")
  subset_table <- data.frame(
    mask = 0:(length(estimated_values$raw) - 1L),
    subset = estimated_values$subset,
    raw = estimated_values$raw,
    projected = projected,
    complement_risk = estimated_values$complement_risk,
    adjusted_risk = estimated_values$adjusted_risk,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      call = call,
      effects = full_effects,
      reduced_effects = reduced_effects,
      selected = selected,
      selected_names = screen$feature_names[selected],
      screen = screen,
      subset_values = subset_table,
      signal = estimated_values$signal,
      noise_floor = estimated_values$noise_floor,
      risk_full = estimated_values$risk_full,
      projection = attr(reduced_effects, "projection"),
      settings = screen$settings
    ),
    class = "prime_fit"
  )
}

#' @export
print.prime_fit <- function(x, ...) {
  retained <- paste(x$selected_names, collapse = ", ")
  if (!nzchar(retained)) retained <- "<none>"
  cat("PRIME fit\n")
  cat("  Retained:", length(x$selected_names), "of", length(x$effects), "\n")
  cat("  Coordinates:", retained, "\n")
  cat("  Use summary() for allocation and diagnostic details.\n")
  invisible(x)
}

#' Summarize a fitted PRIME result
#'
#' @param object A fitted \code{prime_fit} object.
#' @param ... Unused.
#' @return A structured summary of selection, allocation, and diagnostics.
#' @method summary prime_fit
#' @export
summary.prime_fit <- function(object, ...) {
  importance <- sort(
    object$effects[object$selected_names],
    decreasing = TRUE
  )
  importance_table <- matrix(
    unname(importance),
    ncol = 1L,
    dimnames = list(names(importance), "Estimate")
  )
  response_variance <- object$screen$data$response_variance
  explained_fraction <- if (
    !is.null(object$signal) && is.finite(response_variance) &&
      response_variance > 0
  ) {
    object$signal / response_variance
  } else {
    NA_real_
  }
  structure(
    list(
      call = object$call,
      n = object$screen$data$n,
      p = object$screen$data$p,
      candidates = object$screen$candidate_names,
      retained = object$selected_names,
      importance = importance_table,
      allocation_sum = sum(object$effects),
      explained_fraction = explained_fraction,
      signal = object$signal,
      risk_full = object$risk_full,
      projection = object$projection,
      settings = object$settings
    ),
    class = "summary_prime_fit"
  )
}

#' @export
print.summary_prime_fit <- function(x, digits = 4L, ...) {
  retained <- paste(x$retained, collapse = ", ")
  if (!nzchar(retained)) retained <- "<none>"
  cat("Call:\n")
  print(x$call)
  cat("\nData:\n")
  cat(" ", x$n, " observations, ", x$p, " predictors\n", sep = "")
  cat("\nSelection:\n")
  cat(
    " ", length(x$candidates), " candidates, ", length(x$retained),
    " retained (", retained, ")\n",
    sep = ""
  )
  if (nrow(x$importance)) {
    cat("\nPRIME importance:\n")
    print(x$importance, digits = digits)
  }
  cat("\nDiagnostics:\n")
  if (is.finite(x$explained_fraction)) {
    cat(
      " Explained fraction:",
      format(x$explained_fraction, digits = digits),
      "\n"
    )
  }
  if (!is.null(x$projection)) {
    cat(
      " Projection:",
      if (isTRUE(x$projection$converged)) "converged" else "not converged",
      "in", x$projection$iterations, "iteration(s)\n"
    )
  }
  invisible(x)
}

#' Plot PRIME variable importance
#'
#' @param x A fitted \code{prime_fit} object.
#' @param y Unused.
#' @param max_display Maximum number of coordinates to display.
#' @param main Plot title.
#' @param xlab Horizontal-axis label.
#' @param col Color for retained coordinates.
#' @param zero_col Color for excluded coordinates.
#' @param ... Additional graphical parameters passed to \code{plot.default()}.
#' @return The fitted object, invisibly.
#' @method plot prime_fit
#' @export
plot.prime_fit <- function(
  x,
  y = NULL,
  max_display = 20L,
  main = NULL,
  xlab = "PRIME importance",
  col = "#087F8C",
  zero_col = "#A5ADB4",
  ...
) {
  if (!inherits(x, "prime_fit")) {
    stop("x must be a prime_fit object returned by prime().", call. = FALSE)
  }
  if (length(max_display) != 1L || !is.finite(max_display) ||
      max_display < 1 || max_display != as.integer(max_display)) {
    stop("max_display must be a positive integer.", call. = FALSE)
  }

  effects <- x$effects
  order_index <- order(-effects, seq_along(effects))
  order_index <- utils::head(
    order_index,
    min(as.integer(max_display), length(effects))
  )
  values <- effects[order_index]
  labels <- names(values)
  positions <- seq_along(values)
  selected <- labels %in% x$selected_names
  scale <- max(c(values, 0.1))

  graphics::plot.default(
    values,
    positions,
    type = "n",
    xlim = c(-0.025 * scale, 1.12 * scale),
    ylim = c(length(values) + 0.5, 0.5),
    axes = FALSE,
    xlab = xlab,
    ylab = "",
    ...
  )
  ticks <- pretty(c(0, 1.12 * scale))
  graphics::abline(v = ticks, col = "#E4E8EB", lwd = 0.8)
  graphics::segments(0, positions, values, positions, col = "#E4E8EB", lwd = 1.5)
  graphics::points(
    values,
    positions,
    pch = 21,
    bg = ifelse(selected, col, zero_col),
    col = "white",
    cex = ifelse(selected, 1.35, 1.05),
    lwd = 0.8
  )
  graphics::axis(1, at = ticks, col = zero_col, col.axis = "#20262B")
  graphics::axis(2, at = positions, labels = labels, las = 1, tick = FALSE)
  if (!is.null(main)) {
    graphics::title(main = main, adj = 0, line = 0.8)
  }
  invisible(x)
}
