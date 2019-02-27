\name{cv.nclreg_fit}
\alias{cv.nclreg_fit}
\title{Internal function of cross-validation for nclreg}
\description{Internal function to conduct k-fold cross-validation for nclreg, produces a plot,
  and returns cross-validated log-likelihood values for \code{lambda}}
\usage{
cv.nclreg_fit(x, y, weights, lambda=NULL, balance=TRUE, 
rfamily=c("clossR", "closs", "gloss", "qloss"), s=1.5, nfolds=10, foldid, 
type = c("loss", "error"), plot.it=TRUE, se=TRUE, n.cores=2, ...)
}
\arguments{
  \item{x}{\code{x} matrix as in \code{nclreg}.}
  \item{y}{response \code{y} as in \code{nclreg}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
  \item{lambda}{Optional user-supplied lambda sequence; default is
    \code{NULL}, and \code{nclreg} chooses its own sequence}
  \item{balance}{for \code{rfamily="closs", "gloss", "qloss"} only}
  \item{rfamily}{response variable distribution and nonconvex loss function}
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
  \item{\dots}{Other arguments that can be passed to \code{nclreg}.}
}
\details{The function runs \code{nclreg} \code{nfolds}+1 times; the
  first to compute the \code{lambda} sequence, and then to
  compute the fit with each of the folds omitted. The error or the log-likelihood value is
  accumulated, and the average value and standard deviation over the
  folds is computed. Note that \code{cv.nclreg} can be used to search for
  values for \code{alpha}: it is required to call \code{cv.nclreg} with a fixed vector \code{foldid} for different values of \code{alpha}.
}

\value{an object of class \code{"cv.nclreg"} is returned, which is a
  list with the ingredients of the cross-validation fit.
  \item{fit}{a fitted nclreg object for the full data.}
  \item{residmat}{matrix of log-likelihood values with row values for \code{lambda} and column values for \code{k}th cross-validation}
  \item{cv}{The mean cross-validated log-likelihood values - a vector of length
    \code{length(lambda)}.}
  \item{cv.error}{estimate of standard error of \code{cv}.}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in.}
  \item{lambda}{a vector of \code{lambda} values} 
  \item{lambda.which}{index of \code{lambda} that gives minimum \code{cv} value.}
  \item{lambda.optim}{value of \code{lambda} that gives minimum \code{cv} value.}
}

\author{Zhu Wang <wangz1@uthscsa.edu>}
\seealso{\code{\link{nclreg}} and \code{\link{plot}}, \code{\link{predict}}, and \code{\link{coef}} methods for \code{"cv.nclreg"} object.}
\keyword{models}
\keyword{regression}

