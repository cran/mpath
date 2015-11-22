\name{hessianReg}
\alias{hessianReg}
\title{Hessian Matrix of Regularized Estimators}
\description{
 Constructing Hessian matrix for regularized regression parameters.
}
\usage{
hessianReg(x, which, \dots)
}
\arguments{
 \item{x}{a fitted model object.}
 \item{which}{which penalty parameter(s)?}
 \item{\dots}{arguments passed to the \code{meatReg} function.}
}

\details{\code{hessianReg} is a function to compute the Hessian matrix estimate of non-zero regularized estimators. Implemented only for \code{zipath} object with \code{family="negbin"} in the current version. 
}

\value{A matrix containing the Hessian matrix estimate for the non-zero parameters.}

\seealso{\code{\link{breadReg}}, \code{\link{meatReg}}}

\references{
Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}
\author{Zhu Wang <zwang@connecticutchildrens.org>}

\examples{
data("bioChemists", package = "pscl")
fm_zinb <- zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10)
hessianReg(fm_zinb, which=which.min(fm_zinb$bic))
}

\keyword{regression}
