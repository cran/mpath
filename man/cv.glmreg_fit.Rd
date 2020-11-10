\name{cv.glmreg_fit}
\alias{cv.glmreg_fit}
\title{Internal function of cross-validation for glmreg}
\description{Internal function to conduct k-fold cross-validation for glmreg, produces a plot,
  and returns cross-validated log-likelihood values for \code{lambda}}
\usage{
cv.glmreg_fit(x, y, weights, offset, lambda=NULL, balance=TRUE, 
              family=c("gaussian", "binomial", "poisson", "negbin"), 
              type=c("loss", "error"), nfolds=10, foldid, plot.it=TRUE,
              se=TRUE, n.cores=2, trace=FALSE, parallel=FALSE, ...)
}
\arguments{
  \item{x}{\code{x} matrix as in \code{glmreg}.}
  \item{y}{response \code{y} as in \code{glmreg}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of cases. Currently only one offset term can be included in the formula.}
  \item{lambda}{Optional user-supplied lambda sequence; default is
    \code{NULL}, and \code{glmreg} chooses its own sequence}
 \item{balance}{for \code{family="binomial"} only}
 \item{family}{response variable distribution}
\item{type}{ cross-validation criteria. For \code{type="loss"}, loss function (log-negative-likelihood) values and \code{type="error"} is misclassification error if \code{family="binomial"}. } 
 \item{nfolds}{number of folds >=3, default is 10}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in. If supplied,
    \code{nfold} can be missing and will be ignored.}
 \item{plot.it}{ a logical value, to plot the estimated log-likelihood values if \code{TRUE}. }
  \item{se}{ a logical value, to plot with standard errors. }
\item{parallel, n.cores}{a logical value, parallel computing or not with the number of CPU cores to use. The cross-validation loop will attempt to send different CV folds off to different cores.}
 \item{trace}{ a logical value, print progress of cross validation or not} 
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
 \item{fit}{a fitted glmreg object for the full data.}
  \item{residmat}{matrix of log-likelihood values with row values for \code{lambda} and column values for \code{k}th cross-validation}
  \item{cv}{The mean cross-validated log-likelihood values - a vector of
    \code{length(lambda)}.}
 \item{cv.error}{estimate of standard error of \code{cv}.}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in.}
  \item{lambda}{a vector of \code{lambda} values}
  \item{lambda.which}{index of \code{lambda} that gives maximum \code{cv} value.}
  \item{lambda.optim}{value of \code{lambda} that gives maximum \code{cv} value.}
}
\references{
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)
\emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]
}

\author{Zhu Wang <wangz1@uthscsa.edu>}
\seealso{\code{\link{glmreg}} and \code{\link{plot}}, \code{\link{predict}}, and \code{\link{coef}} methods for \code{"cv.glmreg"} object.}
\keyword{models}
\keyword{regression}

