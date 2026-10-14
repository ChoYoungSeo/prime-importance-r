\name{prime_screen}
\alias{prime_screen}
\title{Nonparametric variable selection}
\usage{
prime_screen(X, y, rescale = TRUE, twin_r = NULL, lambda = 1e-6,
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
An object of class \code{prime_screen} with components:
\itemize{
  \item \code{call}: matched function call;
  \item \code{candidates} and \code{candidate_names}: inclusive forward set;
  \item \code{selected} and \code{selected_names}: backward-cleaned set;
  \item \code{forward}: ridge seed, local path, stopping reason, and trace;
  \item \code{backward}: retained set and deletion-comparison trace;
  \item \code{settings}: neighborhood, query, nugget, and calibration rules.
}
}
\description{This function implements the nonparametric model-free variable selection
technique in Cho and Joseph (2026).
}
\details{It constructs a response-relevant candidate set from a BIC-selected prefix of a
standardized ridge ordering and a query-calibrated local-linear forward path, and then
applies sequential query-calibrated backward reduction. Please see Cho and Joseph (2026)
for details.

Search builds an index in each evaluated coordinate subset and retains all
observations in its neighbor pool. If observations are equally distant from a query point, exact search
orders them by their original row numbers. This tie-breaking rule can
choose different neighbors from earlier versions that used the FNN package.
Query contributions are retained for at most two subset fits; allocation caches
scalar risks. Setting \code{neighbor_eps > 0} makes approximation explicit in
\code{settings}. The default preserves exact distances and the neighborhood
size rule.

The ridge ordering is searched through
\eqn{\lceil1.5\log\{\max(2,p)\}\rceil} coordinates, including the intercept-only model,
and BIC selects the initial prefix. The local path adds the coordinate with the
largest explained-variance gain when its paired querywise calibrated gain is
positive. Backward reduction sequentially removes zero calibrated deletion
losses and recomputes all projected neighborhoods after every removal.

Ridge ranking and BIC use a centered unit-variance response internally. The
reported BIC omits the response-scale constant shared by all candidate models;
ridge-path RSS remains in original response units. Paired comparison scales are
computed in response-variance units and returned in original units, avoiding
overflow or underflow of their empirical variances under changes of units.

Use \code{summary(screen)} to see the ridge seed, forward path, retained
coordinates, and backward-reduction counts.
}
\references{
Cho, Y. and Joseph, V. R. (2026). Measuring Variable Importance Using
Proportional Inclusive Marginal Effects. \emph{arXiv:2610.17053}.
\doi{10.48550/arXiv.2610.17053}.
}
\seealso{\code{\link{prime}}}
\examples{
\donttest{
set.seed(5)
X <- matrix(runif(960, -1, 1), nrow = 120)
colnames(X) <- paste0("X", seq_len(ncol(X)))
y <- X[, 1] - X[, 2] + X[, 5]^2 + rnorm(120, sd = 0.08)

screen <- prime_screen(X, y)
screen$selected_names
summary(screen)
screen$forward$trace
screen$backward$trace
}
}
