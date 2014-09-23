\name{cv.glmreg_fit}
\alias{cv.glmreg_fit}
\title{Internal function of cross-validation for glmreg}
\description{Internal function to conduct k-fold cross-validation for glmreg, produces a plot,
  and returns cross-validated loglikelihood values for \code{lambda}}
\usage{
cv.glmreg_fit(x, y, weights, lambda=NULL, balance=TRUE, 
family=c("gaussian", "binomial", "poisson", "negbin"), 
nfolds=10, foldid, plot.it=TRUE, se=TRUE, trace=FALSE, ...)
}
\arguments{
  \item{x}{\code{x} matrix as in \code{glmreg}.}
  \item{y}{response \code{y} as in \code{glmreg}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
  \item{lambda}{Optional user-supplied lambda sequence; default is
    \code{NULL}, and \code{glmreg} chooses its own sequence}
 \item{balance}{for \code{family="binomial"} only}
 \item{family}{response variable distribution}
  \item{nfolds}{number of folds >=3, default is 10}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in. If supplied,
    \code{nfold} can be missing and will be ignored.}
 \item{plot.it}{ a logical value, to plot the estimated loglikelihood values if \code{TRUE}. }
  \item{se}{ a logical value, to plot with standard errors. }
  \item{trace}{ if \code{TRUE}, shows cross-validation progress}
  \item{\dots}{Other arguments that can be passed to \code{glmreg}.}
}
\details{The function runs \code{glmreg} \code{nfolds}+1 times; the
  first to compute the \code{lambda} sequence, and then to
  compute the fit with each of the folds omitted. The error or the log-likelihood value is
  accumulated, and the average value and standard deviation over the
  folds is computed. Note that \code{cv.glmreg} can be used to search for
  values for \code{alpha}: it is required to call \code{cv.glmreg} with a fixed vector \code{foldid} for different values of \code{alpha}.
}

\value{an object of class \code{"cv.glmreg"} is returned, which is a
  list with the ingredients of the cross-validation fit.
  \item{lambda}{the values of \code{lambda} used in the fits.}
  \item{cvm}{The mean cross-validated error - a vector of length
    \code{length(lambda)}.}
  \item{cvsd}{estimate of standard error of \code{cvm}.}
  \item{cvup}{upper curve = \code{cvm+cvsd}.}
  \item{cvlo}{lower curve = \code{cvm-cvsd}.}
  \item{nzero}{number of non-zero coefficients at each \code{lambda}.}
 \item{glmreg.fit}{a fitted glmreg object for the full data.}
  \item{lambda.min}{value of \code{lambda} that gives minimum
    \code{cvm}.}
}
\references{
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)
\emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}, in press
}

\author{Zhu Wang <zwang@connecticutchildrens.org>}
\seealso{\code{\link{glmreg}} and \code{\link{plot}}, \code{\link{predict}}, and \code{\link{coef}} methods for \code{"cv.glmreg"} object.}
\keyword{models}
\keyword{regression}
