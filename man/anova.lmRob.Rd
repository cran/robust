\name{anova.lmRob}
\alias{anova.lmRob}

\title{Compute and Anova table}

\description{
Produces an Anova table for a robustly fitted linear model.
}

\usage{
\method{anova}{lmRob}(object, ..., test = c("RF", "RWald"))
}

\arguments{
\item{object}{an lmRob object.}

\item{\dots}{additional arguments required by the generic anova function.  If \code{\dots} contains additional robustly fitted linear models then the function \code{\link{anova.lmRob.list}} is dispatched.}

\item{test}{a single character value specifying which test should be computed in the Anova table.  The possible choices are "RWald" and "RF".}
}


\details{
The default test used by anova is the \code{"RWald"} test, which is the Wald test based on robust estimates of the coefficients and covariance matrix.  If \code{test} is \code{"RF"}, the robustified F-test is used instead.  
}

\value{
an \code{anova} object.
}

\section{References}{
Hampel, F. R., Ronchetti, E. M., Rousseeuw, P. J., and Stahel, W. A. (1986).  \emph{Robust statistics: the approach based on influence functions}.  John Wiley & Sons.  
}

\seealso{
\code{\link{lmRob}}, 
\code{\link{anova}}, 
\code{\link{anova.lmRob.list}}.
}

\examples{
data(stack.dat)
stack.small <- lmRob(Loss ~ Water.Temp + Acid.Conc., data = stack.dat)
stack.full <- lmRob(Loss ~ ., data = stack.dat)
anova(stack.full)
anova(stack.full, stack.small)
}

\keyword{robust}
\keyword{regression}
\keyword{methods}



