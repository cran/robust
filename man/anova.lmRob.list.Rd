\name{anova.lmRob.list}
\alias{anova.lmRob.list}

\title{Compute and Anova Table}

\description{
Computes an anova table for a list of robustly fitted linear models.  This function is a utility for the user-level \code{\link{anova.lmRob}} function and is not intended to be called directly.
}

\usage{
anova.lmRob.list(object, const, ipsi, yc, test = c("RWald", "RF"), ...)
}

\arguments{
\item{object}{a list of lmRob objects.}

\item{const}{a numeric value containing the tuning constant computed by \code{\link{lmRob.const}}.}

\item{ipsi}{an integer value specifying the psi-function.}

\item{yc}{a numeric value containing the tuning constant computed by \code{\link{lmRob.effvy}}.}

\item{test}{a single character value specifying the test.  Possible values are "RWald" and "RF".}

\item{\dots}{additional arguments required by the generic anova function.}
}

\value{
an \code{anova} object.
}

\seealso{
\code{\link{lmRob}}, 
\code{\link{anova.lmRob}}, 
\code{\link{anova}}.
}

\keyword{methods}
\keyword{robust}


