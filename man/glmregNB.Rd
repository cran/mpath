\name{glmregNB}
\alias{glmregNB}
\alias{glmregNegbin}
\title{fit a negative binomial model with lasso (or elastic net), snet and mnet regularization}
\description{
    Fit a negative binomial linear model via penalized maximum likelihood. The regularization path is computed for the lasso (or elastic net penalty), snet and mnet penalty, at a grid of values for the regularization parameter lambda. }
\usage{
glmregNB(formula, data, weights, offset=NULL, nlambda = 100, lambda=NULL, 
lambda.min.ratio = ifelse(nobs<nvars,0.05,0.001), alpha=1, gamma=3, rescale=TRUE, 
standardize = TRUE, penalty.factor = rep(1, nvars), thresh = 0.001, maxit.theta = 10, 
maxit=1000, eps=.Machine$double.eps, trace=FALSE, start = NULL, 
etastart = NULL, mustart = NULL, theta.fixed=FALSE, theta0=NULL, 
init.theta=ifelse(!theta.fixed, theta0[1],NULL), link=log, 
penalty=c("enet","mnet","snet"), method="glmreg_fit", model=TRUE, 
x.keep=FALSE, y.keep=TRUE, contrasts=NULL, convex=FALSE, ...)
}

\arguments{
  \item{formula}{formula used to describe a model.}
 \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.} 
  \item{weights}{an optional vector of `prior weights' to be used in the fitting process. Should be \code{NULL} or a numeric vector. Default is a vector of 1s with equal weight for each observation.}
   \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor of the model.}
  \item{nlambda}{The number of \code{lambda} values - default is 100.}
 \item{lambda}{A user supplied \code{lambda} sequence}
  \item{lambda.min.ratio}{Smallest value for \code{lambda}, as a fraction of
    \code{lambda.max}, the (data derived) entry value (i.e. the smallest
  value for which all coefficients are zero). The default depends on the
  sample size \code{nobs} relative to the number of variables
  \code{nvars}. If \code{nobs > nvars}, the default is \code{0.001},
  close to zero.  If \code{nobs < nvars}, the default is \code{0.05}.
}
 \item{alpha}{The L2 penalty mixing parameter, with
    \eqn{0\le\alpha\le 1}. \code{alpha=1} is lasso (mcp, scad) penalty; and \code{alpha=0} the ridge penalty.}
\item{gamma}{The tuning parameter of the \code{snet} or \code{mnet} penalty.}
 \item{rescale}{logical value, if TRUE, adaptive rescaling of the penalty parameter for \code{penalty="mnet"} or \code{penalty="snet"} with \code{family} other than "gaussian". See reference}
  \item{standardize}{Logical flag for x variable standardization, prior to
    fitting the model sequence. The coefficients are always returned on
    the original scale. Default is \code{standardize=TRUE}.
  If variables are in the same units already, you might not wish to
  standardize.}
 \item{penalty.factor}{This is a number that multiplies \code{lambda} to allow
    differential shrinkage of coefficients. Can be 0 for some variables, which implies
    no shrinkage, and that variable is always included in the
    model. Default is same shrinkage for all variables.}
  \item{thresh}{Convergence threshold for coordinate descent. Defaults value is \code{1e-6}.}
  \item{maxit.theta}{Maximum number of iterations for estimating \code{theta} scaling parameter}
\item{maxit}{Maximum number of coordinate descent iterations for each \code{lambda} value; default is 1000.}
  \item{eps}{If a number is less than \code{eps} in magnitude, then this number is considered as 0}
  \item{trace}{If \code{TRUE}, fitting progress is reported}
 \item{start, etastart, mustart, \dots}{ arguments for the \code{link{glmreg}} function}
 \item{init.theta}{ initial scaling parameter \code{theta}}
 \item{theta.fixed}{Estimate scale parameter theta? Default is FALSE. Note, the algorithm may become slow. In this case, one may use \code{glmreg} function with \code{family="negbin"}, and a fixed \code{theta}}.
 \item{theta0}{initial scale parameter vector theta, with length \code{nlambda} if \code{theta.fixed=TRUE}. Default is NULL}
\item{convex}{Calculate index for which objective function ceases to
    be locally convex? Default is FALSE and only useful if \code{penalty="mnet" or "snet"}.}
 \item{link}{link function, default is \code{log}}
 \item{penalty}{Type of regularization}
\item{method}{estimation method}
 \item{model, x.keep, y.keep}{logicals. If \code{TRUE} the corresponding components
    of the fit (model frame, response, model matrix) are returned.}
\item{contrasts}{the contrasts corresponding to \code{levels} from the
    respective models}
}
\details{
 The sequence of models implied by \code{lambda} is fit by coordinate
  descent. This is a lasso (mcp, scad) or elastic net (mnet, snet) regularization path
  for fitting the negative binomial linear regression
  paths, by maximizing the penalized log-likelihood.
  Note that the objective function is 
  \deqn{-\sum (weights * loglik) + \lambda*penalty} if \code{standardize=FALSE} and \deqn{-\frac{weights}{\sum(weights)} * loglik + \lambda*penalty} if \code{standardize=TRUE}.
  }
\value{
An object with S3 class \code{"glmreg", "glmregNB"} for the various types of models.
  \item{call}{the call that produced the model fit}
  \item{b0}{Intercept sequence of length \code{length(lambda)}}
  \item{beta}{A \code{nvars x
      length(lambda)} matrix of coefficients.}
  \item{lambda}{The actual sequence of \code{lambda} values used}
  \item{dev}{The computed deviance. The deviance calculations incorporate weights if
  present in the model. The deviance is defined to be 2*(loglike_sat -
  loglike), where loglike_sat is the log-likelihood for the saturated
  model (a model with a free parameter per observation). }
    \item{nulldev}{Null deviance (per observation). This is defined to
  be  2*(loglike_sat -loglike(Null)); The NULL model refers to the
  intercept model.}
  \item{nobs}{number of observations}
}
\references{                                                                          Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)      \emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]
}

\author{Zhu Wang <wangz1@uthscsa.edu>}

\seealso{\code{\link{print}}, \code{\link{predict}}, \code{\link{coef}} and \code{\link{plot}} methods, and the \code{\link{cv.glmregNB}} function.}
\examples{
\dontrun{
data("bioChemists", package = "pscl")
fm_nb <- glmregNB(art ~ ., data = bioChemists)
coef(fm_nb)
### ridge regression
fm <- glmregNB(art ~ ., alpha=0, data = bioChemists, lambda=seq(0.001, 1, by=0.01))
fm <- cv.glmregNB(art ~ ., alpha=0, data = bioChemists, lambda=seq(0.001, 1, by=0.01))
}
}
\keyword{models}
\keyword{regression}
