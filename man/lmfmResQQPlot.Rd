\name{lmfmResQQPlot}
\alias{lmfmResQQPlot}

\title{Residual QQ-Plots with Confidence Envelope}

\description{
Produces side-by-side residual qq-plots from an lmfm object.  A simulated confidence envelope is provided for each plot.
}

\usage{
lmfmResQQPlot(x, envelope = TRUE, half.normal = FALSE, n.samples = 100, level = 0.95, id.n = 3, robustQQline = TRUE, ...)
}

\arguments{
  \item{x}{an lmfm object.}

  \item{envelope}{a logical value.  If TRUE a \code{level} confidence envelope is simulated for each qq-plot.}

  \item{half.normal}{a logical value.  If TRUE the plot is drawn using the absolute values of the residuals.}

  \item{n.samples}{a positive integer value giving the number of samples to compute in the simulation of the confidence envelope.}

  \item{level}{a numeric value between 0 and 1 specifying the confidence level for the envelope.}

  \item{id.n}{a non-negative integer value specifying the number of extreme points to identify.}

  \item{robustQQline}{a logical value.  If TRUE then a qq-line computed using \code{\link{lmRob}} is included in the plot.}

  \item{\dots}{additional arguments are passed to the low-level plotting functions.}
}

\value{
\code{x} is invisibly returned.
}

\keyword{hplot}


