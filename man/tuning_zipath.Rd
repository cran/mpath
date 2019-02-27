\name{tuning.zipath}
\alias{tuning.zipath}

\title{find optimal path for penalized zero-inflated model}

\description{
Fit penalized zero-inflated models, generate multiple paths with varying penalty parameters, therefore determine optimal path with respect to a particular penalty parameter
  }
\usage{
tuning.zipath(formula, data, weights, subset, na.action, offset, standardize=TRUE,
family = c("poisson", "negbin", "geometric"), penalty = c("enet", "mnet", "snet"), 
lambdaCountRatio = .0001, lambdaZeroRatio = c(.1, .01, .001), 
maxit.theta=1, gamma.count=3, gamma.zero=3, ...)

}

\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights. If \code{standardize=TRUE}, weights are renormalized to weights/sum(weights). If \code{standardize=FALSE}, weights are kept as original input}
  \item{subset}{subset of data}
  \item{na.action}{how to deal with missing data}
  \item{offset}{Not implemented yet}
  \item{standardize}{logical value, should variables be standardized?}
  \item{family}{family to fit}
  \item{penalty}{penalty considered as one of \code{enet, mnet, snet}.}
 \item{lambdaCountRatio, lambdaZeroRatio}{Smallest value for \code{lambda.count}
       and \code{lambda.zero}, respectively, as a fraction of
    \code{lambda.max}, the (data derived) entry value (i.e. the smallest
  value for which all coefficients are zero except the intercepts). This lambda.max can be a surrogate value for \code{penalty="mnet"} or \code{"snet"}}
 \item{maxit.theta}{For family="negbin", the maximum iteration allowed for estimating scale parameter theta. Note, the default value 1 is for computing speed purposes, and is typically too small and less desirable in real data analysis}
  \item{gamma.count}{The tuning parameter of the \code{snet} or \code{mnet} penalty for the count part of model.}
  \item{gamma.zero}{The tuning parameter of the \code{snet} or \code{mnet} penalty for the zero part of model.}
  \item{...}{Other arguments passing to \code{zipath}}
}
\details{
From the default \code{lambdaZeroRatio = c(.1, .01, .001)} values, 
find optimal lambdaZeroRatio for penalized zero-inflated Poisson, negative binomial and geometric model 
  }
\value{
  An object of class zipath with the optimal lambdaZeroRatio
}
\references{
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)      \emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]

Zhu Wang, Shuangge Ma, Ching-Yun Wang, Michael Zappitelli, Prasad Devarajan and Chirag R. Parikh (2014)
\emph{EM for Regularized Zero Inflated Regression Models with Applications to Postoperative Morbidity after Cardiac Surgery in Children}, \emph{Statistics in Medicine}. 33(29):5192-208.

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}

\author{Zhu Wang <wangz1@uthscsa.edu>}
\seealso{\code{\link{zipath}}}
\examples{
\dontrun{
## data
data("bioChemists", package = "pscl")

## inflation with regressors
## ("art ~ . | ." is "art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment")
fm_zip2 <- tuning.zipath(art ~ . | ., data = bioChemists, nlambda=10)
summary(fm_zip2)
fm_zinb2 <- tuning.zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10)
summary(fm_zinb2)
}
}


\keyword{models}
\keyword{regression}

 
