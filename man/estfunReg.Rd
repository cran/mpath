\name{estfunReg}
\alias{estfunReg}
\alias{estfunReg.zipath}
\title{Extract Empirical First Derivative of Log-likelihood Function}
\description{
 Generic function for extracting the empirical first derivative of log-likelihood function of a fitted regularized model.
}
\usage{
estfunReg(x, \dots)
}
\arguments{
 \item{x}{a fitted model object.}
 \item{\dots}{arguments passed to methods.}
}

\value{A matrix containing the empirical first derivative of log-likelihood functions.
Typically, this should be an \eqn{n \times k}{n x k} matrix corresponding
to \eqn{n} observations and \eqn{k} parameters. The columns should be named
as in \code{\link{coef}} or \code{\link{terms}}, respectively.
}

\seealso{\code{\link{zipath}}}
\references{
Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}
\author{Zhu Wang <zwang@connecticutchildrens.org>}

\examples{
data("bioChemists", package = "pscl")
fm_zinb <- zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10)
res <- estfunReg(fm_zinb, which=which.min(fm_zinb$bic))
}

\keyword{regression}
