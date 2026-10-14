args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply one installed package library path.")
library_path <- normalizePath(args[1L], mustWork = TRUE)
.libPaths(c(library_path, .libPaths()))
library(prime)
stopifnot(normalizePath(find.package("prime")) == file.path(library_path, "prime"))
set.seed(99113)
n <- 20000L
X <- matrix(rnorm(n * 6L), nrow=n)
y <- X[, 1L]^2 + rnorm(n)
queries <- rep.int(seq_len(n) - 1L, 10L)
cat("READY_FOR_OWN_PROCESS_SIGINT\n")
flush(stdout())
tryCatch({
  prime:::indexed_forward_values_cpp(
    X, y, queries, c(0L, 1L), 2:5, 213L,
    lambda=1e-6, response_variance=var(y), threads=2L,
    eps=0, query_block_size=16L
  )
  cat("UNEXPECTED_COMPLETION\n")
  quit(save="no", status=3L)
}, interrupt=function(e) {
  cat("CAUGHT_INTERRUPT_CLEANLY\n")
  quit(save="no", status=0L)
}, error=function(e) {
  cat("UNEXPECTED_ERROR:", conditionMessage(e), "\n")
  quit(save="no", status=2L)
})
