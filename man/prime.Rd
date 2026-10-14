\name{prime}
\alias{prime}
\title{Proportional Inclusive Marginal Effects}
\usage{
prime(X, y, rescale = TRUE, twin_r = NULL, lambda = 1e-6,
  neighbor_eps = 0, query_block_size = 256L)
}
\arguments{
  \item{X}{Numeric predictor matrix or data frame with unique column names.}
  \item{y}{Finite numeric response vector with length \code{nrow(X)}.}
  \item{rescale}{Logical; standardize predictors before distance calculations.}
  \item{twin_r}{Optional integer inverse Data-Twinning split ratio. \code{NULL}
  uses every observation as a query location.}
  \item{lambda}{Nonnegative slope-only numerical nugget used in local-linear
  fits.}
  \item{neighbor_eps}{Nonnegative relative Euclidean distance tolerance.
  Zero uses exact indexed search. Positive values enable approximate neighbors;
  the distance tolerance does not bound allocation or selection error.}
  \item{query_block_size}{Positive query batch size. The native engine streams
  individual neighborhoods and checks for interrupts at these intervals.}
}
\value{
An object of class \code{prime_fit} with components:
\itemize{
  \item \code{call}: matched function call;
  \item \code{effects}: full-length, zero-padded PRIME allocation;
  \item \code{reduced_effects}: allocation on retained coordinates;
  \item \code{selected} and \code{selected_names}: retained coordinates;
  \item \code{screen}: complete \code{prime_screen} result;
  \item \code{subset_values}: raw and projected subset values with complement
  and noise-adjusted risks;
  \item \code{signal}, \code{noise_floor}, and \code{risk_full}: estimated
  explained variation and full retained-set risk;
  \item \code{projection}: monotone-projection diagnostics;
  \item \code{settings}: fixed numerical and statistical settings.
}
If screening retains no coordinate, \code{effects} is a zero vector and
\code{subset_values} is \code{NULL}, indicating that a normalized allocation
is unavailable.
}
\description{
Estimates Proportional Inclusive Marginal Effects of Cho and Joseph (2026).
}
\details{The function executes BIC-gated ridge-start screening, query-calibrated forward and
backward reduction, and proportional allocation (Herin et al. 2024) on the retained total-effect
set function. Please see Cho and Joseph (2026) for details.

Exact allocation requires \eqn{2^d} subset values and \eqn{d2^{d-1}}
projection edges for \eqn{d} retained coordinates. The edge-index capacity is
checked before subset evaluation; dimensions 28 and above are rejected by
the current 32-bit projection indices. Lower dimensions can still be costly
in memory and runtime.

Search builds an index in each evaluated coordinate subset and retains all
observations in its neighborhood pool. If observations are equally distant from a query point, exact search
orders them by their original row numbers. This tie-breaking rule can
choose different neighbors from earlier versions that used the FNN package.
Query contributions are retained for at most two subset fits; allocation caches
scalar risks. Setting \code{neighbor_eps > 0} makes approximation explicit in
\code{settings}. The default preserves exact distances and the neighborhood
size rule.

For retained dimension \eqn{d}, the local neighborhood size is
\deqn{K_d=\left\lceil\max\{1.5\sqrt n,\,3d\}\right\rceil,}
subject to the local-fit feasibility cap. The local nugget is applied to slopes
only. Every observation is used as a query by default; Data Twinning reduces
the outer query average without removing observations from the neighbor pool.

Exact allocation evaluates \eqn{2^d} subset values and therefore scales
exponentially with the retained dimension.

Use \code{summary(fit)} to see retained coordinates, importance estimates,
explained variation, and projection diagnostics.
}
\references{
Herin, M., Il Idrissi, M., Chabridon, V., and Iooss, B. (2024).
Proportional marginal effects for global sensitivity analysis.
\emph{SIAM/ASA Journal on Uncertainty Quantification}, 12(2), 667--692.

Cho, Y. and Joseph, V. R. (2026). Measuring Variable Importance Using
Proportional Inclusive Marginal Effects. \emph{arXiv:2610.17053}.
\doi{10.48550/arXiv.2610.17053}.
}
\seealso{\code{\link{plot.prime_fit}}, \code{\link{prime_screen}}}
\examples{
\donttest{
set.seed(5)
X <- matrix(runif(960, -1, 1), nrow = 120)
colnames(X) <- paste0("X", seq_len(ncol(X)))
y <- X[, 1] - X[, 2] + X[, 5]^2 + rnorm(120, sd = 0.08)

fit <- prime(X, y)
fit$selected_names
round(fit$effects, 3)
summary(fit)
plot(fit)
}
}
