\name{sandwichReg}
\alias{sandwichReg}
\title{Making Sandwiches with Bread and Meat for Regularized Estimators}
\description{
 Constructing sandwich covariance matrix estimators by
 multiplying bread and meat matrices for regularized regression parameters.
}
\usage{
sandwichReg(x, breadreg.=breadReg, meatreg.=meatReg, which, log=FALSE, \dots)
}
\arguments{
 \item{x}{a fitted model object.}
 \item{breadreg.}{either a breadReg matrix or a function for computing
   this via \code{breadreg.(x)}.}
 \item{meatreg.}{either a breadReg matrix or a function for computing
   this via \code{meatreg.(x, \dots)}.}
 \item{which}{which penalty parameters(s) to compute?}
 \item{log}{if TRUE, the corresponding element is with respect to log(theta) in negative binomial regression. Otherwise, for theta} 
 \item{\dots}{arguments passed to the \code{meatReg} function.}
}

\details{\code{sandwichReg} is a function to compute an estimator for the covariance of the non-zero parameters. It takes a breadReg matrix (i.e., estimator of the expectation of the negative
derivative of the penalized estimating functions) and a meatReg matrix (i.e.,
estimator of the variance of the log-likelihood function) and multiplies
them to a sandwich with meat between two slices of bread. By default
\code{\link{breadReg}} and \code{\link{meatReg}} are called. Implemented only for \code{zipath} object with \code{family="negbin"} in the current version. 
}

\value{A matrix containing the sandwich covariance matrix estimate for the non-zero parameters.}

\seealso{\code{\link{breadReg}}, \code{\link{meatReg}}}

\references{
Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}
\author{Zhu Wang <wangz1@uthscsa.edu>}

\examples{
data("bioChemists", package = "pscl")
fm_zinb <- zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10, maxit.em=1)
sandwichReg(fm_zinb, which=which.min(fm_zinb$bic))
}

\keyword{regression}
