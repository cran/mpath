\name{cv.glmreg}
\alias{cv.glmreg}
\alias{cv.glmreg.formula}
\alias{cv.glmreg.matrix}
\alias{cv.glmreg.default}
\alias{plot.cv.glmreg}
\alias{coef.cv.glmreg}
\title{Cross-validation for glmreg}
\description{Does k-fold cross-validation for glmreg, produces a plot,
  and returns cross-validated loglikelihood values for \code{lambda}}
\usage{
\method{cv.glmreg}{formula}(formula, data, weights, offset=NULL, ...)
\method{cv.glmreg}{matrix}(x, y, weights, offset=NULL, ...)
\method{cv.glmreg}{default}(x,  ...)
\method{plot}{cv.glmreg}(x,se=TRUE,ylab=NULL, main=NULL, width=0.02, col="darkgrey", ...)
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
\references{                                                                                              Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)      \emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}, in press
}
\author{Zhu Wang <zwang@connecticutchildrens.org>}

\seealso{\code{\link{glmreg}} and \code{\link{plot}}, \code{\link{predict}}, and \code{\link{coef}} methods for \code{"cv.glmreg"} object.}
\examples{
\dontrun{
data("bioChemists", package = "pscl")
fm_pois <- cv.glmreg(art ~ ., data = bioChemists, family = "poisson")
plot(fm_pois)
title("Poisson Family",line=2.5)
}
}
\keyword{models}
\keyword{regression}

