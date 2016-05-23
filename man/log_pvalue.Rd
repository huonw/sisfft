\name{log_pvalue}
\alias{log_pvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Logarithm of a tail sum of the \code{L}-fold convolution of
  \code{exp(log_pmf)} with bounded relative error.
}
\description{
  \code{log_pvalue} returns the logarithm of an approximation to the
  tail sum from index \code{s0} of the
  convolution of \code{exp(log_pmf)} with itself \code{L} times. The
  approximation has relative error bounded by \code{gamma}. That
  is, if \code{P = sum(q[s0:length(q)])} with \code{q = p * p * ... *
    p} (where the \code{*} denote \code{L} convolutions), then the
  approximation \code{Q} satisfies \code{abs(Q - P) <= gamma * P}.
}
\usage{
log_pvalue(log_pmf, s0, L, gamma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{log_pmf}{
    a numeric vector, the logarithm of a pmf.
}
  \item{s0}{
    a positive integer, the index for which the tail sum should be computed.
  }
  \item{L}{
    a positive integer giving the the number of times to convolve
    \code{log_pmf} with itself.
}
  \item{gamma}{
    a positive number less than 0.5, to control the relative error.
  }
}
\value{
  A number that is the logarithm of the approximation to the tail sum.
}
\references{
%% ~put references to the literature/web site here ~
}

\examples{
pmf <- c(0.5, 0.5)
log_pval <- log_pvalue(log(pmf), 52, 101, 1e-3)
abs(exp(log_pval) - 0.5) <= 1e-3 * 0.5

pmf <- exp(30:1) / sum(exp(30:1))
# 2901 is the last entry in the convolution
log_pval <- log_pvalue(log(pmf), 2901, 100, 1e-3)
# (much) smaller than smallest double
log_pval < log(.Machine$double.xmin)

exact <- log(pmf[30]) * 100
# relative error, without underflow
abs(expm1(log_pval - exact)) <= 1e-3
}
\seealso{
  \code{\link{log_convolve}}
  \code{\link{log_convolve_power}}
}

\keyword{htest}
