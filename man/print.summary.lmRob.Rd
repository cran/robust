\name{print.summary.lmRob}
\alias{print.summary.lmRob}

\title{Print Method for summary.lmRob Class Objects}

\description{
Displays the formula, coefficients, degrees of freedom, and the residual standard error in a fitted lmRob object.
}

\usage{
print.summary.lmRob(x, ...)
}

\arguments{
  \item{x}{a summary.lmRob object.}

  \item{\dots}{additional arguments required by the generic print function.  In particular the 'digits' argument is useful for specifying the number of significant digits when displaying numeric output.}
}

\value{
\code{x} is invisibly returned.
}

\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{
\code{\link{lmRob}}, 
\code{\link{lmRob.object}}.
}

\examples{
data(stack.dat)
stack.robsum <- summary(lmRob(Loss ~ ., data = stack.dat))
print(stack.robsum, digits = 4)
}

\keyword{methods}


