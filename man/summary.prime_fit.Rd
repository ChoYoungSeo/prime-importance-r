\name{summary.prime_fit}
\alias{summary.prime_fit}
\alias{print.summary_prime_fit}
\keyword{internal}
\title{Summarize a Fitted PRIME Result}
\usage{
\method{summary}{prime_fit}(object, ...)

\method{print}{summary_prime_fit}(x, digits = 4L, ...)
}
\arguments{
  \item{object}{A fitted \code{prime_fit} object returned by \code{prime()}.}
  \item{x}{A \code{summary_prime_fit} object.}
  \item{digits}{Number of significant digits used when printing.}
  \item{...}{Additional arguments for the generic methods.}
}
\value{
\code{summary.prime_fit()} returns a structured list containing the matched
call, sample and predictor dimensions, candidate and retained coordinates, a
PRIME importance estimate matrix, explained fraction, allocation total,
projection diagnostics, and settings. The print method returns this object
invisibly.
}
\description{
Builds and prints a structured summary of a complete PRIME fit.
}
\seealso{\code{\link{prime}}, \code{\link{plot.prime_fit}}}
