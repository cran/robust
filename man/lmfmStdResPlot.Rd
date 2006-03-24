\name{lmfmStdResPlot}
\alias{lmfmStdResPlot}

\title{Standardized Residual Plot}

\description{
Produces side-by-side standardized residual plots from an lmfm object.
}

\usage{
lmfmStdResPlot(x, level = 0.95, id.n = 3, ...)
}

\arguments{
  \item{x}{an lmfm object.}

  \item{level}{a numeric value between 0 and 1 specifying the confidence level used to draw the threshold in the plot.}

  \item{id.n}{a non-negative integer value specifying the number of extreme points to identify.}

  \item{\dots}{additional arguments are passed to the low-level plotting functions.}
}

\value{
\code{x} is invisibly returned.
}

\keyword{hplot}


