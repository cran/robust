\name{lmfmSRvsRDPlot}
\alias{lmfmSRvsRDPlot}

\title{Standardized Residuals versus Robust Distances Plot}

\description{
Produces side-by-side residuals versus robust distances plots from an lmfm object.
}

\usage{
lmfmSRvsRDPlot(x, level = 0.95, id.n = 3, ...)
}

\arguments{
  \item{x}{an lmfm object.}

  \item{level}{a numeric value between 0 and 1 giving the confidence level used to draw the threshold for both the distances and residuals axes.}

  \item{id.n}{a non-negative integer value specifying the number of extreme points to identify.}

  \item{\dots}{additional arguments are passed to the low-level plotting functions.}
}

\value{
\code{x} is invisibly returned.
}

\keyword{hplot}


