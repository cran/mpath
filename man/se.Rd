\name{se}
\alias{se}
\alias{se.zipath}
\title{Standard Error of Regularized Estimators}
\description{
 Generic function for computing standard errors of non-zero regularized estimators
}
\usage{
se(x, which, log=TRUE, \dots)
}
\arguments{
 \item{x}{a fitted model object.}
 \item{which}{which penalty parameter(s)?}
 \item{log}{if TRUE, the computed standard error is for log(theta) for negative binomial regression, otherwise, for theta.}
 \item{\dots}{arguments passed to methods.}
}

\value{A vector containing standard errors of non-zero regularized estimators.}

\seealso{\code{\link{zipath}}}
\references{
Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}
\author{Zhu Wang <wangz1@uthscsa.edu>}

\examples{
data("bioChemists", package = "pscl")
fm_zinb <- zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10, maxit.em=1)
res <- se(fm_zinb, which=which.min(fm_zinb$bic))
}

\keyword{regression}
