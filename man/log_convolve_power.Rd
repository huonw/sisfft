\name{log_convolve_power}
\alias{log_convolve_power}

\title{
  Logarithm of the \code{L}-fold convolution of \code{exp(log_pmf)},
  with bounded relative error.
}
\description{
  \code{log_convolve_power} returns the logarithm of an approximation
 of the convolution of \code{exp(log_pmf)} with itself \code{L}
 times. The approximation has relative error bounded by \code{alpha},
  along with the subtractive factor of \code{delta}. If \code{p} is the exact convolution and \code{q} is the
 approximation, then it satisfies \code{p[i] * (1 - 1 / alpha) -
 delta <= q[i] <= p[i] * (1 + 1 / alpha)} for all \code{i}.
}

\usage{
log_convolve_power(log_pmf, L, alpha, delta)
}

\arguments{
  \item{log_pmf}{a numeric vector, the logarithm of a pmf.}
  \item{L}{a positive integer giving the the number of times to convolve
    \code{log_pmf} with itself.}
  \item{alpha}{a number larger than 2, to control the relative error.}
  \item{delta}{a non-negative number, to act as a lower bound on
  elements needed}
}
\value{A vector that is the logarithm of the approximation to the convolution.}
\examples{
 pmf <- c(0.5, 0.5)
 convolve <- log_convolve_power(log(pmf), 3, 1e3, 0)
 exact <- c(0.125, 0.375, 0.375, 0.125)
 (exact * (1 - 1e-3) <= pmf) & (pmf <= exact * (1 + 1e-3))
}

\seealso{
  \code{\link{log_convolve}}
  \code{\link{log_pvalue}}
}

\keyword{array}
