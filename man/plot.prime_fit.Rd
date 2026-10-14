\name{plot.prime_fit}
\alias{plot.prime_fit}
\title{Plot PRIME Variable Importance}
\usage{
\method{plot}{prime_fit}(x, y = NULL, max_display = 20L,
  main = NULL, xlab = "PRIME importance",
  col = "#087F8C", zero_col = "#A5ADB4", ...)
}
\arguments{
  \item{x}{A fitted \code{prime_fit} object returned by \code{prime()}, such
  as \code{obj} in the example below.}
  \item{y}{Optional argument retained for compatibility with the generic
  \code{plot(x, y, ...)} interface; no second input is needed.}
  \item{max_display}{Maximum number of coordinates to display.}
  \item{main}{Optional plot title.}
  \item{xlab}{Horizontal-axis label.}
  \item{col}{Color for retained coordinates.}
  \item{zero_col}{Color for excluded coordinates.}
  \item{...}{Additional graphical parameters passed to \code{plot.default()}.}
}
\value{
The fitted object, invisibly.
}
\description{
Displays PRIME importance from a fitted result. Retained
coordinates are highlighted and excluded coordinates remain visible at zero.
}
\seealso{\code{\link{prime}}}
\examples{
\donttest{
set.seed(5)
X <- matrix(runif(960, -1, 1), nrow = 120)
colnames(X) <- paste0("X", seq_len(ncol(X)))
y <- X[, 1] - X[, 2] + X[, 5]^2 + rnorm(120, sd = 0.08)

obj <- prime(X, y)
plot(obj)
}
}
