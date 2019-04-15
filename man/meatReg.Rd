\name{meatReg}
\alias{meatReg}
\title{Meat Matrix Estimator}
\description{
 Estimating the variance of the first derivative of log-likelihood function}
 
\usage{
meatReg(x, which, \dots)
}
\arguments{
 \item{x}{a fitted model object. Currently only implemented for \code{zipath} object with \code{family="negbin"}}
 \item{which}{which penalty parameter(s)?}
 \item{\dots}{arguments passed to the \code{\link{estfunReg}} function.}
}

\details{
See reference below
}

\value{A \deqn{k \times k}{k x k} covariance matrix of first derivative of log-likelihood function}

\seealso{\code{\link{sandwichReg}}, \code{\link{breadReg}}, \code{\link{estfunReg}}}

\references{
Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}
\author{Zhu Wang <wangz1@uthscsa.edu>}

\examples{
data("bioChemists", package = "pscl")
fm_zinb <- zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10, maxit.em=1)
meatReg(fm_zinb, which=which.min(fm_zinb$bic))
}

\keyword{regression}
