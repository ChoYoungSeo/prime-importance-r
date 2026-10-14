args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply one installed package library path.")
library_path <- normalizePath(args[1L], mustWork = TRUE)
.libPaths(c(library_path, .libPaths()))
library(prime)
stopifnot(normalizePath(find.package("prime")) == file.path(library_path, "prime"))
options(warn = 2)
x <- seq(-2, 2, length.out = 60L)
X <- cbind(x, x, 2*x, sin(3*x), cos(2*x))
y <- x^2 + sin(7*x)/5
queries <- c(0L, 10L, 30L, 59L, 10L)
selected <- c(1L, 2L)
candidates <- 3:5
k <- 15L
expected <- vapply(candidates, function(candidate) {
  Xs <- X[, sort(c(selected, candidate)), drop = FALSE]
  sse <- vapply(queries + 1L, function(q) {
    d2 <- rowSums(sweep(Xs, 2L, Xs[q, ], "-")^2)
    eligible <- setdiff(seq_len(nrow(Xs)), q)
    neighbors <- eligible[order(d2[eligible], eligible)][seq_len(k)]
    design <- cbind(1, sweep(Xs[neighbors, , drop = FALSE], 2L,
                            Xs[q, ], "-"))
    decomp <- svd(design)
    keep <- decomp$d > max(dim(design)) * .Machine$double.eps * max(decomp$d)
    U <- decomp$u[, keep, drop = FALSE]
    sum((y[neighbors] - U %*% crossprod(U, y[neighbors]))^2)
  }, numeric(1))
  (var(y) - mean(sse)/(k - ncol(Xs) - 1L))/var(y)
}, numeric(1))
cat("prime", as.character(packageVersion("prime")), "\n")
for (workers in c(1L, 2L)) {
  actual <- prime:::indexed_forward_values_cpp(
    X, y, queries, selected-1L, candidates-1L, k,
    lambda=0, response_variance=var(y), threads=workers,
    eps=0, query_block_size=2L
  )$values
  max_error <- max(abs(as.numeric(actual)-expected))
  stopifnot(all(is.finite(actual)), max_error < 1e-10)
  cat("lambda=0 workers=", workers, " max absolute score error=",
      format(max_error, digits=16), "\n", sep="")
}
