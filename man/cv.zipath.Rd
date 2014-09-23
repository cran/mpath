\name{cv.zipath}
\alias{cv.zipath}
\alias{coef.cv.zipath}
\title{Cross-validation for zipath}
\description{Does k-fold cross-validation for zipath, produces a plot,
  and returns cross-validated loglikelihood values for \code{lambda}}
\usage{
cv.zipath(formula, data, weights, nlambda=100, lambda.count=NULL, lambda.zero=NULL,
nfolds=10, foldid, plot.it=TRUE, se=TRUE, trace=FALSE,...)
\method{coef}{cv.zipath}(object, which=object$lambda.which, model = c("full", "count", "zero"), \dots)
}
\arguments{
  \item{formula}{symbolic description of the model}
  \item{data}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
\item{nlambda}{number of \code{lambda} value, default value is 10.}
  \item{lambda.count}{Optional user-supplied lambda.count sequence; default is
    \code{NULL}}
  \item{lambda.zero}{Optional user-supplied lambda.zero sequence; default is
    \code{NULL}}
  \item{nfolds}{number of folds >=3, default is 10}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in. If supplied,
    \code{nfold} can be missing and will be ignored.}
 \item{plot.it}{ a logical value, to plot the estimated loglikelihood values if \code{TRUE}. }
  \item{se}{ a logical value, to plot with standard errors. }
  \item{trace}{ if \code{TRUE}, shows cross-validation progress }
  \item{\dots}{Other arguments that can be passed to \code{zipath}.}
 \item{object}{object of class \code{cv.zipath}.}
 \item{which}{Indices of the pair of penalty parameters \code{lambda.count} and \code{lambda.zero} at which
    estimates are extracted.  By default, the one which generates the optimal cross-validation value.}
 \item{model}{character specifying for which component of the model the
    estimated coefficients should be extracted.}
}
\details{The function runs \code{zipath} \code{nfolds}+1 times; the
  first to compute the \code{(lambda.count, lambda.zero)} sequence, and then to
  compute the fit with each of the folds omitted. The log-likelihood value is
  accumulated, and the average value and standard deviation over the
  folds is computed. Note that \code{cv.zipath} can be used to search for
  values for \code{count.alpha} or \code{zero.alpha}: it is required to call \code{cv.zipath} with a fixed vector \code{foldid} for different values of \code{count.alpha} or \code{zero.alpha}.

The method for \code{\link[stats]{coef}} by default
  return a single vector of coefficients, i.e., all coefficients are concatenated. By setting the \code{model}
  argument, the estimates for the corresponding model components can be extracted.
}
\value{an object of class \code{"cv.zipath"} is returned, which is a
  list with the components of the cross-validation fit.
 \item{fit}{a fitted zipath object for the full data.}
 \item{residmat}{matrix for cross-validated log-likelihood at each \code{(count.lambda, zero.lambda)} sequence}
  \item{bic}{The cross-validated BIC value - a vector of length
    \code{length(count.lambda)}.}
 \item{fraction}{a sequence from 1:\code{nlambda}. \code{nlambda} is the same as the argument if any one of \code{(count.lambda, zero.lambda)} is missing; otherwise \code{nlambda=length(count.lambda)}}
  \item{cv}{The mean cross-validated log-likelihood - a vector of length
    \code{length(count.lambda)}.}
  \item{cv.error}{estimate of standard error of \code{cv}.}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in.} 
 \item{lambda.which}{index of \code{(count.lambda, zero.lambda)} that gives minimum \code{cv}.}
 \item{lambda.optim}{value of \code{(count.lambda, zero.lambda)} that gives minimum \code{cv}.}
}
\references{
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)      \emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}, in press

Zhu Wang, Shuangge Ma, Ching-Yun Wang, Michael Zappitelli, Prasad Devarajan and Chirag R. Parikh (2014)
\emph{EM for Regularized Zero Inflated Regression Models with Applications to Postoperative Morbidity after Cardiac Surgery in Children}, \emph{Statistics in Medicine}, in press

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2014) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, manuscript
}

\author{Zhu Wang <zwang@connecticutchildrens.org>}

\seealso{\code{\link{zipath}} and \code{\link{plot}}, \code{\link{predict}}, and \code{\link{coef}} methods for \code{"cv.zipath"} object.}
\examples{
\dontrun{
data("bioChemists", package = "pscl")
fm_zip <- cv.zipath(art ~ . | ., data = bioChemists, family = "poisson", nlambda=10)
coef(fm_zip)
fm_znb <- cv.zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10)
coef(fm_znb)
}
}
\keyword{models}
\keyword{regression}
