\name{summary.glmregNB}
\alias{summary.glmregNB}
\alias{print.summary.glmregNB}
\title{
Summary Method Function for Objects of Class 'glmregNB'
}
\description{
Summary results of fitted penalized negative binomial regression model
}
\usage{
\method{summary}{glmregNB}(object, \dots)
}
\arguments{
\item{object}{
fitted model object of class \code{glmregNB}.
}
\item{\dots}{
arguments passed to or from other methods.
}}
\value{
Summary of fitted penalized negative binomial model
}
\details{
This function is a method for the generic function
\code{summary()} for class \code{"glmregNB"}.
It can be invoked by calling \code{summary(x)} for an
object \code{x} of the appropriate class, or directly by
calling \code{summary.glmregNB(x)} regardless of the
class of the object.
}
\references{
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)
\emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]
}

\author{Zhu Wang <zwang@connecticutchildrens.org>}
\seealso{
\code{\link{summary}}, \code{\link[MASS]{glm.nb}}}
\keyword{models}
\examples{
\dontrun{
data(quine, package="MASS")
summary(glmregNB(Days ~ Eth*Age*Lrn*Sex, quine, link = log))
}
}
