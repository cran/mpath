\name{cv.glmreg}
\alias{cv.glmreg}
\alias{cv.glmreg.formula}
\alias{cv.glmreg.matrix}
\alias{cv.glmreg.default}
\alias{plot.cv.glmreg}
\alias{predict.cv.glmreg}
\alias{coef.cv.glmreg}
\title{Cross-validation for glmreg}
\description{Does k-fold cross-validation for glmreg, produces a plot,
  and returns cross-validated log-likelihood values for \code{lambda}}
\usage{
\method{cv.glmreg}{formula}(formula, data, weights, offset=NULL, ...)
\method{cv.glmreg}{matrix}(x, y, weights, offset=NULL, ...)
\method{cv.glmreg}{default}(x,  ...)
\method{plot}{cv.glmreg}(x,se=TRUE,ylab=NULL, main=NULL, width=0.02, col="darkgrey", ...)
\method{predict}{cv.glmreg}(object, newx, ...)
\method{coef}{cv.glmreg}(object,which=object$lambda.which, ...)
}
\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{x}{\code{x} matrix as in \code{glmreg}. It could be object of \code{cv.glmreg}.}
  \item{y}{response \code{y} as in \code{glmreg}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
  \item{offset}{Not implemented yet}
  \item{object}{object of \code{cv.glmreg}}
  \item{newx}{Matrix of values at which predictions are to be made.  Not
    used for \code{type="coefficients"}}
  \item{which}{Indices of the penalty parameter \code{lambda} at which
    estimates are extracted. By default, the one which generates the optimal cross-validation value.}
  \item{se}{logical value, if \code{TRUE}, standard error curve is also plotted} 
  \item{ylab}{ylab on y-axis}
  \item{main}{title of plot}
  \item{width}{width of lines}
  \item{col}{color of standard error curve}
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
  \item{bic}{matrix of BIC values with row values for \code{lambda} and column values for \code{k}th cross-validation}
  \item{cv}{The mean cross-validated log-likelihood values - a vector of length
    \code{length(lambda)}.}
 \item{cv.error}{estimate of standard error of \code{cv}.}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in.}
 \item{lambda}{a vector of \code{lambda} values} 
  \item{lambda.which}{index of \code{lambda} that gives maximum \code{cv} value.}
  \item{lambda.optim}{value of \code{lambda} that gives maximum \code{cv} value.}
}
\references{Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)      \emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]
}
\author{Zhu Wang <wangz1@uthscsa.edu>}

\seealso{\code{\link{glmreg}} and \code{\link{plot}}, \code{\link{predict}}, and \code{\link{coef}} methods for \code{"cv.glmreg"} object.}
\examples{
data("bioChemists", package = "pscl")
fm_pois <- cv.glmreg(art ~ ., data = bioChemists, family = "poisson")
title("Poisson Family",line=2.5)
predict(fm_pois, newx=bioChemists)
coef(fm_pois)
}
\keyword{models}
\keyword{regression}

