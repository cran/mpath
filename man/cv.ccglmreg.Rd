\name{cv.ccglmreg}
\alias{cv.ccglmreg}
\alias{cv.ccglmreg.formula}
\alias{cv.ccglmreg.matrix}
\alias{cv.ccglmreg.default}
\alias{plot.cv.ccglmreg}
\alias{coef.cv.ccglmreg}
\title{Cross-validation for ccglmreg}
\description{Does k-fold cross-validation for ccglmreg, produces a plot,
  and returns cross-validated log-likelihood values for \code{lambda}}
\usage{
\method{cv.ccglmreg}{formula}(formula, data, weights, offset=NULL, ...)
\method{cv.ccglmreg}{matrix}(x, y, weights, offset=NULL, ...)
\method{cv.ccglmreg}{default}(x,  ...)
\method{plot}{cv.ccglmreg}(x,se=TRUE,ylab=NULL, main=NULL, width=0.02, col="darkgrey", ...)
\method{coef}{cv.ccglmreg}(object,which=object$lambda.which, ...)
}
\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{x}{\code{x} matrix as in \code{ccglmreg}. It could be object of \code{cv.ccglmreg}.}
  \item{y}{response \code{y} as in \code{ccglmreg}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
  \item{offset}{Not implemented yet}
  \item{object}{object of \code{cv.ccglmreg}}
  \item{which}{Indices of the penalty parameter \code{lambda} at which
    estimates are extracted. By default, the one which generates the optimal cross-validation value.}
  \item{se}{logical value, if \code{TRUE}, standard error curve is also plotted} 
  \item{ylab}{ylab on y-axis}
  \item{main}{title of plot}
  \item{width}{width of lines}
  \item{col}{color of standard error curve}
  \item{\dots}{Other arguments that can be passed to \code{ccglmreg}.}
}
\details{The function runs \code{ccglmreg} \code{nfolds}+1 times; the
  first to compute the \code{lambda} sequence, and then to
  compute the fit with each of the folds omitted. The error or the loss value is
  accumulated, and the average value and standard deviation over the
  folds is computed. Note that \code{cv.ccglmreg} can be used to search for
  values for \code{alpha}: it is required to call \code{cv.ccglmreg} with a fixed vector \code{foldid} for different values of \code{alpha}. 
}
\value{an object of class \code{"cv.ccglmreg"} is returned, which is a
  list with the ingredients of the cross-validation fit.
  \item{fit}{a fitted ccglmreg object for the full data.}
  \item{residmat}{matrix of log-likelihood values with row values for \code{lambda} and column values for \code{k}th cross-validation}
  \item{bic}{matrix of BIC values with row values for \code{lambda} and column values for \code{k}th cross-validation}
  \item{cv}{The mean cross-validated log-likelihood values - a vector of length
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

