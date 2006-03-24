\name{lmfmResVsFittedPlot}
\alias{lmfmResVsFittedPlot}

\title{Residuals versus Fitted Values Plot}

\description{
Produces side-by-side residuals versus fitted values plots from an lmfm object.
}

\usage{
lmfmResVsFittedPlot(x, smooths = FALSE, rugplot = FALSE, id.n = 3, ...)
}

\arguments{
  \item{x}{an lmfm object.}

  \item{smooths}{a logical value.  If TRUE then a smooth curve computed using \code{\link{loess}} is included in each panel of the plot.}

  \item{rugplot}{a logical value.  If TRUE then a rugplot is included in each panel of the plot.}

  \item{id.n}{a non-negative integer value giving the number of extreme points to identify.}

  \item{\dots}{additional arguments are passed to the low-level plotting functions.}
}

\value{
\code{x} is invisibly returned.
}

\keyword{utilities}


