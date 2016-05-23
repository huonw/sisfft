\name{log_convolve}
\alias{log_convolve}
\title{
  Logarithm of the convolution of \code{exp(log_pmf1)} and \code{exp(log_pmf2)},
  with bounded relative error.
}
\description{
  \code{log_convolve} returns the logarithm of an approximation of the
 convolution of \code{pmf1 = exp(log_pmf1)} with \code{pmf2 = exp(log_pmf2)}. The
 approximation has relative error bounded by \code{beta}, after
 accounting for the lower bound \code{delta}. Denote \code{p} the exact
 convolution of \code{pmf1} and \code{pmf2}, \code{q} the exact
 convolution of those vectors after setting any value less than
 \code{delta} to zero, and finally \code{v} the approximation. The
 approximation satisfies
 \code{(1 - beta) * q[i] <= v[i] <= (1 + beta) * q[i]} for all \code{i}.
}
\usage{
log_convolve(log_pmf1, log_pmf2, beta, delta = 0)
}
\arguments{
  \item{log_pmf1}{
    a vector that is the logarithm of a pmf.
}
  \item{log_pmf2}{
    a vector that is the logarithm of a pmf.
}
  \item{beta}{
    a number less than 0.5, to control the relative error. Note that
    this differs to Wilson, Keich 2016 which describes aFFT-C with
    \code{alpha := 1 / beta}.
}
  \item{delta}{
    a non-negative number, to act as a lower bound on elements needed.
}
}
\value{A vector that is the logarithm of the approximation to the convolution.}
\references{
Huon Wilson, Uri Keich, Accurate pairwise convolutions of non-negative
vectors via FFT, Computational Statistics & Data Analysis, Volume 101,
September 2016, Pages 300-315, ISSN 0167-9473, http://dx.doi.org/10.1016/j.csda.2016.03.010.
}

\seealso{
  \code{\link{log_convolve_power}}
  \code{\link{log_pvalue}}
}
\examples{
pmf1 <- c(0.5, 0.5)
pmf2 <- c(0.1, 0.2, 0.3, 0.4)

log_conv <- log_convolve(log(pmf1), log(pmf2), 1e-3, 0)
conv <- exp(log_conv)

exact <- c(0.05, 0.15, 0.25, 0.35, 0.20)
(exact * (1 - 1e-3) < conv) & (conv < exact * (1 + 1e-3))

# the same convolution with a bound that allows dropping the 0.1
log_conv_bound <- log_convolve(log(pmf1), log(pmf2), 1e-3, 0.15)
conv_bound <- exp(log_conv_bound)

exact_bound <- c(0.00, 0.10, 0.25, 0.35, 0.20)
(exact_bound * (1 - 1e-3) < conv_bound) & (conv_bound < exact * (1 + 1e-3))
}
\keyword{array}
