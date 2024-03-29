\name{cv.zipath_fit}
\alias{cv.zipath_fit}
\title{Cross-validation for zipath}
\description{Internal function k-fold cross-validation for zipath, produces a plot, and returns cross-validated log-likelihood values for \code{lambda}}
\usage{
cv.zipath_fit(X, Z, Y, weights, offsetx, offsetz, nlambda=100, lambda.count=NULL,
              lambda.zero=NULL, nfolds=10, foldid, plot.it=TRUE, se=TRUE, 
              n.cores=2, trace=FALSE, parallel=FALSE, ...)
}
\arguments{
   \item{X}{predictor matrix of the count model}
   \item{Z}{predictor matrix of the zero model}
   \item{Y}{response variable}
   \item{weights}{optional numeric vector of weights.}
   \item{offsetx}{optional numeric vector with an a priori known component to be included in the linear predictor of the count model.}
   \item{offsetz}{optional numeric vector with an a priori known component to be included in the linear predictor of the zero model.}
  \item{nlambda}{number of \code{lambda} value, default value is 10.}
  \item{lambda.count}{Optional user-supplied lambda.count sequence; default is
    \code{NULL}}
  \item{lambda.zero}{Optional user-supplied lambda.zero sequence; default is
    \code{NULL}}
  \item{nfolds}{number of folds >=3, default is 10}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in. If supplied,
    \code{nfold} can be missing and will be ignored.}
 \item{plot.it}{ a logical value, to plot the estimated log-likelihood values if \code{TRUE}. }
  \item{se}{ a logical value, to plot with standard errors. }
 \item{n.cores}{The number of CPU cores to use. The cross-validation loop
  will attempt to send different CV folds off to different cores.}
\item{trace}{ a logical value, print progress of cross-validation or not }
 \item{parallel}{ a logical value, parallel computing or not}
  \item{\dots}{Other arguments that can be passed to \code{zipath}.}
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
  \item{bic}{matrix of BIC values with row values for \code{lambda} and column values for \code{k}th cross-validation}
  \item{cv}{The mean cross-validated log-likelihood - a vector of length
    \code{length(count.lambda)}.}
  \item{cv.error}{estimate of standard error of \code{cv}.}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in.} 
 \item{lambda.which}{index of \code{(count.lambda, zero.lambda)} that gives maximum \code{cv}.}
 \item{lambda.optim}{value of \code{(count.lambda, zero.lambda)} that gives maximum \code{cv}.}
}
\references{
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)      \emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]

Zhu Wang, Shuangge Ma, Ching-Yun Wang, Michael Zappitelli, Prasad Devarajan and Chirag R. Parikh (2014)
\emph{EM for Regularized Zero Inflated Regression Models with Applications to Postoperative Morbidity after Cardiac Surgery in Children}, \emph{Statistics in Medicine}. 33(29):5192-208.

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}

\author{Zhu Wang <zwang145@uthsc.edu>}

\seealso{\code{\link{zipath}} and \code{\link{plot}}, \code{\link{predict}}, and \code{\link{coef}} methods for \code{"cv.zipath"} object.}
\keyword{models}
\keyword{regression}

