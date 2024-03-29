\name{lognormRob}
\alias{lognormRob}

\title{Robust Estimation of Lognormal Distribution Parameters}

\description{
  Robust estimation of lognormal distribution parameters.
}

\usage{
lognormRob(x, estim = c("tdmean"), control = lognormRob.control(estim, \dots), \dots)
}
\arguments{
  \item{x}{a numeric vector containing the sample.}
  \item{estim}{a character string specifying which estimator to use.}
  \item{control}{a list of control parameters appropriate for the estimator in \code{estim}.}
  \item{\dots}{control parameters may also be given here.}
}

\value{
  a list with class \dQuote{fitdstn} containing the following elements:

  \item{estimate}{a named numeric vector containing the parameter estimates.}
  \item{sd}{a named numeric vector containing the standard deviations of
    the parameter estimates.  Missing in current implementation.}
  \item{vcov}{a numeric matrix containing the variance-covariance matrix
    of the estimated parameter vector.  Missing in current
    implementation.}
  \item{mu}{a single numeric value containing an estimate of the mean.}
  \item{V.mu}{a single numeric value containing the variance of the estimated mean.}
  \item{control}{a list containing the control parameters used by the estimator.}
}
\note{
  The \code{print} method displays the estimated parameters and their
  standard errors (in parentheses).
}

\seealso{
  \code{\link{lognormRob.control}}, \code{\link{fitdstnRob}}.
}

\keyword{robust}


