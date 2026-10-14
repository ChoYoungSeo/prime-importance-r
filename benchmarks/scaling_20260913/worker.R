#!/usr/bin/env Rscript
# One package, one dataset, and one operation per fresh R process.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(11L, 12L))
out <- normalizePath(args[1L])
writeLines(as.character(Sys.getpid()), file.path(out, "worker.pid"))
library_path <- normalizePath(args[2L])
variant <- args[3L]
mode <- args[4L]
n <- as.integer(args[5L]); p <- as.integer(args[6L]); seed <- as.integer(args[7L])
neighbor_eps <- as.numeric(args[8L]); query_block_size <- as.integer(args[9L])
max_selected <- as.integer(args[10L]); recall_queries <- as.integer(args[11L])
twin_r <- if (length(args) == 11L || args[12L] == "NULL") NULL else as.integer(args[12L])
if (!is.null(twin_r)) stopifnot(is.finite(twin_r), twin_r >= 2L, twin_r <= floor(n / 2L))
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(prime))
stopifnot(normalizePath(find.package("prime")) == file.path(library_path, "prime"))
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
writeLines(.libPaths(), file.path(out, "library_paths.txt"))
writeLines(capture.output(Sys.getenv(c(
  "PRIME_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
  "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "R_LIBS_USER", "R_LIBS"
))), file.path(out, "thread_environment.txt"))
internal <- function(name) getFromNamespace(name, "prime")
stage <- function(label) writeLines(label, file.path(out, "stage.txt"))
stage("generate data")
# The response and active coordinates do not change when nuisance p changes.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
X <- vapply(seq_len(p), function(j) {
  set.seed(seed + 1000L * j)
  runif(n, -1, 1)
}, numeric(n))
colnames(X) <- paste0("X", seq_len(p))
set.seed(seed + 900000L)
y <- 2 * sin(pi * X[, 1L]) + 1.5 * X[, 2L] + X[, 3L]^2 + rnorm(n, sd = 0.2)
saveRDS(list(X = X, y = y), file.path(out, "input.rds"), compress = FALSE)
input_md5 <- unname(tools::md5sum(file.path(out, "input.rds")))
# Fingerprint the exact inputs, then remove the large duplicate artifact.
unlink(file.path(out, "input.rds"))
metrics <- list(status = "completed", error = "", input_md5 = input_md5,
                actual_package_version = as.character(packageVersion("prime")),
                r_version = R.version.string, n = n, p = p,
                query_count = n, twin_r = if (is.null(twin_r)) "NULL" else twin_r,
                neighbor_eps = neighbor_eps,
                screen_sec = NA_real_, pme_sec = NA_real_, operation_sec = NA_real_)
extra <- if (variant == "baseline") list() else list(
  neighbor_eps = neighbor_eps, query_block_size = query_block_size
)
screen <- NULL
fit <- NULL
started <- proc.time()[["elapsed"]]
tryCatch({
  if (mode == "neighbors") {
    stage("neighbor recall")
    q <- unique(as.integer(round(seq(1, n, length.out = min(n, recall_queries)))))
    k <- min(n - 1L, as.integer(ceiling(1.5 * sqrt(n))))
    # Limit metric dimension and query count; this diagnostic returns q by k.
    Xn <- scale(X[, seq_len(min(5L, p)), drop = FALSE])
    begin <- proc.time()[["elapsed"]]
    indexed <- internal("indexed_neighbors_cpp")(Xn, q - 1L, k, neighbor_eps)
    metrics$neighbor_sec <- proc.time()[["elapsed"]] - begin
    ids <- as.matrix(indexed$indices) + 1L
    squared <- as.matrix(indexed$squared_distances)
    stopifnot(identical(dim(ids), c(length(q), k)), all(ids >= 1L & ids <= n))
    recall <- ratio <- numeric(length(q))
    distance_error <- 0
    for (i in seq_along(q)) {
      distance <- rowSums((sweep(Xn, 2L, Xn[q[i], ], "-"))^2)
      distance[q[i]] <- Inf
      reference <- head(order(distance, seq_len(n)), k)
      stopifnot(!q[i] %in% ids[i, ], !anyDuplicated(ids[i, ]))
      recall[i] <- length(intersect(reference, ids[i, ])) / k
      ratio[i] <- sqrt(max(distance[ids[i, ]]) / distance[reference[k]])
      distance_error <- max(distance_error, max(abs(distance[ids[i, ]] - squared[i, ])))
    }
    metrics$neighbor_recall_mean <- mean(recall)
    metrics$neighbor_recall_min <- min(recall)
    metrics$neighbor_distance_ratio_max <- max(ratio)
    metrics$neighbor_squared_distance_error_max <- distance_error
    metrics$neighbor_queries <- length(q)
    metrics$neighbor_k <- k
    write.csv(data.frame(query = q, recall = recall, kth_distance_ratio = ratio),
              file.path(out, "neighbor_recall.csv"), row.names = FALSE)
  } else {
    if (mode %in% c("screen", "pipeline")) {
      stage("screen")
      begin <- proc.time()[["elapsed"]]
      screen <- do.call(prime::prime_screen, c(list(X = X, y = y, twin_r = twin_r), extra))
      metrics$screen_sec <- proc.time()[["elapsed"]] - begin
      selected <- screen$selected
      metrics$candidates <- paste(screen$candidates, collapse = ",")
      metrics$candidate_size <- length(screen$candidates)
      metrics$forward_steps <- nrow(screen$forward$trace)
      metrics$backward_steps <- nrow(screen$backward$trace)
      metrics$forward_stop <- screen$forward$stop_reason
      saveRDS(list(forward = screen$forward, backward = screen$backward,
                   settings = screen$settings), file.path(out, "screen_diagnostics.rds"))
      write.csv(screen$forward$trace, file.path(out, "forward.csv"), row.names = FALSE)
      write.csv(screen$backward$trace, file.path(out, "backward.csv"), row.names = FALSE)
      data <- screen$data
    } else {
      # Isolate the PME estimator at the same fixed s=3 in every variant.
      stage("prepare fixed active-set PME")
      selected <- 1:3
      data <- do.call(internal(".prime_prepare"), c(list(
        X = X, y = y, rescale = TRUE, twin_r = twin_r, lambda = 1e-6
      ), extra))
      metrics$candidates <- "1,2,3"
      metrics$candidate_size <- 3L
    }
    metrics$query_count <- length(data$query_index)
    writeLines(as.character(data$query_index), file.path(out, "query_indices.txt"))
    metrics$selected <- paste(selected, collapse = ",")
    metrics$selected_size <- length(selected)
    metrics$active_inclusion <- all(1:3 %in% selected)
    metrics$false_positives <- sum(!selected %in% 1:3)
    if (mode %in% c("pme", "pipeline")) {
      if (length(selected) > max_selected) {
        metrics$status <- "allocation_dimension_budget_exceeded"
      } else {
        stage("PME subset estimation and allocation")
        begin <- proc.time()[["elapsed"]]
        effects <- numeric(p)
        if (length(selected)) {
          values <- internal(".prime_estimate_subset_values")(data, selected)
          allocation <- prime::prime_allocate(values$raw)
          effects[selected] <- allocation
          metrics$signal <- values$signal
          metrics$noise_floor <- values$noise_floor
          metrics$coalitions <- length(values$raw)
          metrics$projection_converged <- attr(allocation, "projection")$converged
          write.csv(data.frame(mask = seq_along(values$raw) - 1L, raw = values$raw,
                              projected = attr(allocation, "projected_values")),
                    file.path(out, "coalitions.csv"), row.names = FALSE)
        }
        metrics$pme_sec <- proc.time()[["elapsed"]] - begin
        metrics$allocation_sum <- sum(effects)
        write.csv(data.frame(coordinate = seq_len(p), effect = effects),
                  file.path(out, "effects.csv"), row.names = FALSE)
      }
    }
  }
}, error = function(condition) {
  metrics$status <<- "error"
  metrics$error <<- conditionMessage(condition)
})
metrics$operation_sec <- proc.time()[["elapsed"]] - started
write.csv(as.data.frame(metrics, stringsAsFactors = FALSE),
          file.path(out, "metrics.csv"), row.names = FALSE)
stage(metrics$status)
if (metrics$status == "error") quit(status = 1L)
