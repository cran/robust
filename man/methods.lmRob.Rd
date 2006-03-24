\name{methods.lmRob}
\alias{coef.lmRob}
\alias{deviance.lmRob}
\alias{formula.lmRob}
\alias{labels.lmRob}
\alias{model.frame.lmRob}
\alias{model.matrix.lmRob}
\alias{residuals.lmRob}
\alias{weights.lmRob}

\title{
Methods for lmRob Generic Functions
}

\description{
Accessor methods for the coefficients, deviance, formula, labels, model frame, model matrix, residuals, and weights in an \code{\link{lmRob}} object.
}

\usage{
coef.lmRob(object, ...)
deviance.lmRob(object, ...)
formula.lmRob(x, ...)
labels.lmRob(object, ...)
model.frame.lmRob(formula, ...)
model.matrix.lmRob(object, ...)
residuals.lmRob(object, ...)
weights.lmRob(object, ...)
}

\arguments{
\item{formula}{an lmRob object.}

\item{object}{an lmRob object.}

\item{x}{an lmRob object.}

\item{...}{additional arguments required by the generic functions.}
}

\value{
A vector, matrix or formula containing the requested component of the lmRob object.
}

\seealso{
\code{\link{lmRob}}, 
\code{\link{lmRob.object}}.
}

\examples{
data(stack.dat)
stack.rob <- lmRob(Loss ~ ., data = stack.dat) 
coef.lmRob(stack.rob) 
deviance.lmRob(stack.rob) 
formula.lmRob(stack.rob) 
labels.lmRob(stack.rob) 
model.frame.lmRob(stack.rob) 
model.matrix.lmRob(stack.rob) 
residuals.lmRob(stack.rob) 
weights.lmRob(stack.rob) 
}

\keyword{methods}


