\name{cv.ccglmreg_fit}
\alias{cv.ccglmreg_fit}
\title{Internal function of cross-validation for ccglmreg}
\description{Internal function to conduct k-fold cross-validation for ccglmreg, produces a plot,
  and returns cross-validated loss values for \code{lambda}}
\usage{
cv.ccglmreg_fit(x, y, weights, offset, lambda=NULL, balance=TRUE, cfun=4, dfun=1, 
                s=1.5, nfolds=10, foldid, type = c("loss", "error"), plot.it=TRUE,
                se=TRUE, n.cores=2, trace=FALSE, parallel=FALSE, ...)
}
\arguments{
  \item{x}{\code{x} matrix as in \code{ccglmreg}.}
  \item{y}{response \code{y} as in \code{ccglmreg}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of cases. Currently only one offset term can be included in the formula.}
  \item{lambda}{Optional user-supplied lambda sequence; default is
    \code{NULL}, and \code{ccglmreg} chooses its own sequence}
  \item{balance}{for \code{dfun=4, 5, 6} only}
  \item{cfun}{ a number from 1 to 7, type of convex cap (concave) function }
  \item{dfun}{ a number from 1, 4-7, type of convex downward function}
  \item{s}{ nonconvex loss tuning parameter for robust regression and classification. }
  \item{nfolds}{number of folds >=3, default is 10}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in. If supplied,
    \code{nfold} can be missing and will be ignored.}
  \item{type}{ cross-validation criteria. For \code{type="loss"}, loss function values and \code{type="error"} is misclassification error. }
  \item{plot.it}{ a logical value, to plot the estimated log-likelihood values if \code{TRUE}. }
  \item{se}{ a logical value, to plot with standard errors. }
  \item{n.cores}{The number of CPU cores to use. The cross-validation loop
    will attempt to send different CV folds off to different cores.}
  \item{trace}{ a logical value, print progress of cross validation or not}
  \item{parallel}{ a logical value, parallel computing or not}
  \item{\dots}{Other arguments that can be passed to \code{ccglmreg}.}
}
\details{The function runs \code{ccglmreg} \code{nfolds}+1 times; the
  first to compute the \code{lambda} sequence, and then to
  compute the fit with each of the folds omitted. The error or the log-likelihood value is
  accumulated, and the average value and standard deviation over the
  folds is computed. Note that \code{cv.ccglmreg} can be used to search for
  values for \code{alpha}: it is required to call \code{cv.ccglmreg} with a fixed vector \code{foldid} for different values of \code{alpha}.
}

\value{an object of class \code{"cv.ccglmreg"} is returned, which is a
  list with the ingredients of the cross-validation fit.
  \item{fit}{a fitted ccglmreg object for the full data.}
  \item{residmat}{matrix of loss values or errors with row values for \code{lambda} and column values for \code{k}th cross-validation}
  \item{cv}{The mean cross-validated loss values or errors - a vector of length
    \code{length(lambda)}.}
  \item{cv.error}{estimate of standard error of \code{cv}.}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in.}
  \item{lambda}{a vector of \code{lambda} values} 
  \item{lambda.which}{index of \code{lambda} that gives minimum \code{cv} value.}
  \item{lambda.optim}{value of \code{lambda} that gives minimum \code{cv} value.}
}
\references{
   Zhu Wang (2020)
   \emph{Unified Robust Estimation via the COCO, arXiv e-prints}, \url{https://arxiv.org/abs/2010.02848}
  }

\author{Zhu Wang <wangz1@uthscsa.edu>}
\seealso{\code{\link{ccglmreg}} and \code{\link{plot}}, \code{\link{predict}}, and \code{\link{coef}} methods for \code{"cv.ccglmreg"} object.}
\keyword{models}
\keyword{regression}

