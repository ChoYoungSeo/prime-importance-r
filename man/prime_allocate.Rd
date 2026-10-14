\name{prime_allocate}
\alias{prime_allocate}
\keyword{internal}
\title{Allocate Precomputed Subset Values}
\usage{
prime_allocate(subset_values, zero_tolerance = sqrt(.Machine$double.eps))
}
\arguments{
  \item{subset_values}{Finite numeric vector ordered by binary subset masks from
  zero through \eqn{2^d-1}. Its length must be \eqn{2^d}.}
  \item{zero_tolerance}{Finite nonnegative numeric threshold relative to the
  full-set value, applied to proper nonempty subsets after projection.
  Zero retains positive projected values.}
}
\value{
A named numeric allocation. The importance values sum to the projected
full-set value. Attributes contain:
\itemize{
  \item \code{projected_values}: the nonnegative monotone projection after
  applying the relative null threshold, in the supplied game units;
  \item \code{projection}: iteration count, convergence flag, and maximum
  remaining violation.
}
}
\description{
Projects supplied subset values onto the set of
nonnegative monotone set functions and evaluates the canonical
proportional-value extension by ratio-potential recursion.
}
\details{
The first value is set to zero. The full-set value must be nonnegative and is
preserved. Positive-full games are normalized internally, then their allocations
and projected values are returned in the supplied units. Ratio potentials are
evaluated in log space. A zero-full game projects to zero and returns zero
allocations. Input ratios that overflow or underflow during normalization are
rejected with a numerical-range error.

The projection convergence tolerance does not snap positive values to zero.
The relative null threshold affects only proper nonempty subsets; the empty
and full endpoints stay fixed. Thus a threshold of one or greater produces
the unanimity game when the full-set value is positive. Projection diagnostics
report violations in the supplied game units.

Computation and storage scale exponentially with dimension
\eqn{d}.
Subset masks permit at most 30 coordinates. The current exact projection also
requires its \eqn{d2^{d-1}} Hasse edges to fit 32-bit integer indices, which
rejects dimensions 28 and above before projection. Passing this representation
check does not guarantee practical memory use or runtime.
}
\references{
Herin, M., Il Idrissi, M., Chabridon, V., and Iooss, B. (2024).
Proportional marginal effects for global sensitivity analysis.
\emph{SIAM/ASA Journal on Uncertainty Quantification}, 12(2), 667--692.
\doi{10.1137/22M153032X}.
}
\examples{
values <- c(0, 0.25, 0.60, 1)
allocation <- prime_allocate(values)
allocation
attr(allocation, "projected_values")
}
