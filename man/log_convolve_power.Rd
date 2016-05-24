\name{log_convolve_power}
\alias{log_convolve_power}

\title{
  Logarithm of the \code{L}-fold convolution of \code{exp(log_pmf)},
  with bounded relative error.
}
\description{
  \code{log_convolve_power} returns the logarithm of an approximation
 of the convolution of \code{exp(log_pmf)} with itself \code{L}
 times. The approximation has relative error bounded by \code{beta},
  along with the subtractive factor of \code{delta}. If \code{p} is the exact convolution and \code{v} is the
 approximation, then it satisfies \code{p[i] * (1 - beta) -
 delta <= v[i] <= p[i] * (1 + beta)} for all \code{i}.
}

\usage{
log_convolve_power(log_pmf, L, beta, delta)
}

\arguments{
  \item{log_pmf}{a numeric vector, the logarithm of a pmf.}
  \item{L}{a positive integer giving the the number of times to convolve
    \code{log_pmf} with itself.}
  \item{beta}{a number less than 0.5, to control the relative error.}
  \item{delta}{a non-negative number, to act as a lower bound on
  elements needed.}
}
\value{A vector that is the logarithm of the approximation to the
  convolution.}

\references{
Huon Wilson, Uri Keich, Accurate small tail probabilities of sums of iid
lattice-valued random variables via FFT, Submitted.
}

\examples{
 pmf <- c(0.5, 0.5)
 convolve <- log_convolve_power(log(pmf), 3, 1e-3, 0)
 exact <- c(0.125, 0.375, 0.375, 0.125)
 (exact * (1 - 1e-3) <= pmf) & (pmf <= exact * (1 + 1e-3))
}

\seealso{
  \code{\link{log_convolve}}
  \code{\link{log_pvalue}}
}

\keyword{array}
