\name{summary.prime_screen}
\alias{summary.prime_screen}
\alias{print.summary_prime_screen}
\keyword{internal}
\title{Summarize PRIME Screening}
\usage{
\method{summary}{prime_screen}(object, ...)

\method{print}{summary_prime_screen}(x, ...)
}
\arguments{
  \item{object}{A \code{prime_screen} object returned by
  \code{prime_screen()}.}
  \item{x}{A \code{summary_prime_screen} object.}
  \item{...}{Additional arguments for the generic methods.}
}
\value{
\code{summary.prime_screen()} returns a structured list containing the matched
call, sample and predictor dimensions, ridge seed, local path, candidate and
retained coordinates, stopping reason, backward evaluation and deletion
counts, and settings. The print method returns this object invisibly.
}
\description{
Builds and prints a structured summary of the PRIME forward and backward
selection stages.
}
\seealso{\code{\link{prime_screen}}, \code{\link{summary.prime_fit}}}
