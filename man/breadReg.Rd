\name{breadReg}
\alias{breadReg}
\alias{breadReg.zipath}
\title{Bread for Sandwiches in Regularized Estimators}
\description{
 Generic function for extracting an estimator for the bread of
 sandwiches.
}
\usage{
breadReg(x, which, \dots)
}
\arguments{
 \item{x}{a fitted model object.}
 \item{which}{which penalty parameter(s)?}
 \item{\dots}{arguments passed to methods.}
}

\value{A matrix containing an estimator for the penalized second derivative of log-likelihood function.
Typically, this should be an \eqn{k \times k}{k x k} matrix corresponding
to \eqn{k} parameters. The rows and columns should be named
as in \code{\link{coef}} or \code{\link{terms}}, respectively.}

\seealso{\code{\link{meatReg}}, \code{\link{sandwichReg}}}

\references{

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.

}
\author{Zhu Wang <zwang@connecticutchildrens.org>}

\examples{
data("bioChemists", package = "pscl")
fm_zinb <- zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10)
breadReg(fm_zinb, which=which.min(fm_zinb$bic))
}
\keyword{regression}
