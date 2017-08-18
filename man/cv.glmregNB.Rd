\name{cv.glmregNB}
\alias{cv.glmregNB}
\title{Cross-validation for glmregNB}
\description{Does k-fold cross-validation for glmregNB, produces a plot,
  and returns cross-validated log-likelihood values for \code{lambda}}
\usage{
cv.glmregNB(formula, data, weights, lambda=NULL,
nfolds=10, foldid, plot.it=TRUE, se=TRUE, n.cores=2, ...)
}
\arguments{
  \item{formula}{symbolic description of the model}
  \item{data}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
  \item{lambda}{Optional user-supplied lambda sequence; default is
    \code{NULL}, and \code{glmregNB} chooses its own sequence}
  \item{nfolds}{number of folds - default is 10. Although \code{nfolds}
    can be as large as the sample size (leave-one-out CV), it is not
    recommended for large datasets. Smallest value allowable is \code{nfolds=3}}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in. If supplied,
    \code{nfold} can be missing.}
 \item{plot.it}{ a logical value, to plot the estimated log-likelihood values if \code{TRUE}. }
  \item{se}{ a logical value, to plot with standard errors. }
 \item{n.cores}{The number of CPU cores to use. The cross-validation loop
  will attempt to send different CV folds off to different cores.}
  \item{\dots}{Other arguments that can be passed to \code{glmregNB}.}
}
\details{The function runs \code{glmregNB} \code{nfolds}+1 times; the
  first to get the \code{lambda} sequence, and then the remainder to
  compute the fit with each of the folds omitted. The error is
  accumulated, and the average error and standard deviation over the
  folds is computed. 
  Note that \code{cv.glmregNB} does NOT search for
  values for \code{alpha}. A specific value should be supplied, else
  \code{alpha=1} is assumed by default. If users would like to
  cross-validate \code{alpha} as well, they should call \code{cv.glmregNB}
  with a pre-computed vector \code{foldid}, and then use this same fold vector
  in separate calls to \code{cv.glmregNB} with different values of
  \code{alpha}. 
}
\value{an object of class \code{"cv.glmregNB"} is returned, which is a
  list with the ingredients of the cross-validation fit.
 \item{fit}{a fitted glmregNB object for the full data.}
  \item{residmat}{matrix of log-likelihood values with row values for \code{lambda} and column values for \code{k}th cross-validation}
  \item{cv}{The mean cross-validated log-likelihood values - a vector of length
    \code{length(lambda)}.}
  \item{cv.error}{The standard error of cross-validated log-likelihood values - a vector of length
    \code{length(lambda)}.}
\item{lambda}{a vector of \code{lambda} values} 
  \item{foldid}{indicators of data used in each cross-validation, for reproductive purposes}
  \item{lambda.which}{index of \code{lambda} that gives maximum \code{cv} value.}
  \item{lambda.optim}{value of \code{lambda} that gives maximum \code{cv} value.}
}
\references{
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)
\emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print] 
}
\author{Zhu Wang <zwang@connecticutchildrens.org>}
\seealso{\code{\link{glmregNB}} and \code{\link{plot}}, \code{\link{predict}}, and \code{\link{coef}} methods for \code{"cv.glmregNB"} object.}
\examples{
\dontrun{
data("bioChemists", package = "pscl")
fm_nb <- cv.glmregNB(art ~ ., data = bioChemists)
plot(fm_nb)
}
}
\keyword{models}
\keyword{regression}

