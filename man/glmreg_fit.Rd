\name{glmreg_fit}
\alias{glmreg_fit}
\title{Internal function to fit a GLM with lasso (or elastic net), snet and mnet regularization}
\description{
    Fit a generalized linear model via penalized maximum likelihood.  The
  regularization path is computed for the lasso (or elastic net penalty), snet and mnet penalty, at a grid
  of values for the regularization parameter lambda. Fits
  linear, logistic, Poisson and negative binomial (fixed scale parameter) regression models.
  }
\usage{
glmreg_fit(x, y, weights, start=NULL, etastart=NULL, mustart=NULL, 
nlambda=100, lambda=NULL, lambda.min.ratio=ifelse(nobs<nvars,.05, .001),alpha=1, 
gamma=3, rescale=TRUE, standardize=TRUE, penalty.factor = rep(1, nvars),thresh=1e-6, 
eps.bino=1e-5, maxit=1000, eps=.Machine$double.eps, theta, 
family=c("gaussian", "binomial", "poisson", "negbin"), penalty=c("enet","mnet","snet"), 
convex=FALSE, x.keep=FALSE, y.keep=TRUE, trace=FALSE)
}

\arguments{
  \item{x}{input matrix, of dimension nobs x nvars; each row is an
  observation vector.} 
  \item{y}{response variable. Quantitative for \code{family="gaussian"}.
  Non-negative counts for \code{family="poisson"} or \code{family="negbin"}. For
  \code{family="binomial"} should be either a factor with two levels or
  a vector of proportions.} 
  \item{weights}{observation weights. Can be total counts if responses are proportion matrices. Default is 1 for each observation}
  \item{start}{starting values for the parameters in the linear predictor.}
  \item{etastart}{starting values for the linear predictor.}
  \item{mustart}{starting values for the vector of means.}
  \item{nlambda}{The number of \code{lambda} values - default is 100. The sequence may be truncated before \code{nlambda} is reached if a close to saturated model is fitted. See also \code{satu}.}
  \item{lambda}{ by default, the algorithm provides a sequence of regularization values, or a user supplied \code{lambda} sequence }
  \item{lambda.min.ratio}{Smallest value for \code{lambda}, as a fraction of
    \code{lambda.max}, the (data derived) entry value (i.e. the smallest
  value for which all coefficients are zero except the intercept). Note, there is no closed formula for \code{lambda.max} in general. If \code{rescale=TRUE}, \code{lambda.max} is the same for \code{penalty="mnet" or "snet"}. Otherwise, some modifications are required. For instance, for small \code{gamma} value, half of the square root (if \code{lambda.max} is too small) of the computed \code{lambda.max} can be used when \code{penalty="mnet" or "snet"}. The default of \code{lambda.min.ratio} depends on the
  sample size \code{nobs} relative to the number of variables
  \code{nvars}. If \code{nobs > nvars}, the default is \code{0.001},
  close to zero.  If \code{nobs < nvars}, the default is \code{0.05}.}
  \item{alpha}{The \eqn{L_2} penalty mixing parameter, with
    \eqn{0<\le alpha\le 1}. \code{alpha=1} is lasso (mcp, scad) penalty; and \code{alpha=0} the ridge penalty. However, if \code{alpha=0}, one must provide \code{lambda} values.}
 \item{gamma}{The tuning parameter of the \code{snet} or \code{mnet} penalty.}
 \item{rescale}{logical value, if TRUE, adaptive rescaling of the penalty parameter for \code{penalty="mnet"} or \code{penalty="snet"} with \code{family} other than "gaussian". See reference}
  \item{standardize}{logical value for x variable standardization, prior to
    fitting the model sequence. The coefficients are always returned on
    the original scale. Default is \code{standardize=TRUE}.}
  \item{penalty.factor}{This is a number that multiplies \code{lambda} to allow
    differential shrinkage of coefficients. Can be 0 for some variables, which implies
    no shrinkage, and that variable is always included in the
    model. Default is same shrinkage for all variables.}
  \item{thresh}{Convergence threshold for coordinate descent. Defaults value is \code{1e-6}.}
  \item{eps.bino}{a lower bound of probabilities to be claimed as zero, for computing weights and related values when \code{family="binomial"}.}
    \item{maxit}{Maximum number of coordinate descent iterations for each \code{lambda} value; default is 1000.}
    \item{eps}{If a coefficient is less than \code{eps} in magnitude, then it is reported to be 0}
 \item{convex}{Calculate index for which objective function ceases to
    be locally convex? Default is FALSE and only useful if \code{penalty="mnet" or "snet"}.}
 \item{theta}{ an overdispersion scaling parameter for \code{family="negbin"}}
  \item{family}{Response type (see above)}
 \item{penalty}{Type of regularization}
 \item{x.keep, y.keep}{ For glmreg: logical values indicating whether the response
          vector and model matrix used in the fitting process should be
          returned as components of the returned value.
          For glmreg_fit: x is a design matrix of dimension n * p,
          and x is a vector of observations of length n.}
 \item{trace}{If \code{TRUE}, fitting progress is reported}
}
\details{
  The sequence of models implied by \code{lambda} is fit by coordinate
  descent. For \code{family="gaussian"} this is the lasso, mcp or scad sequence if
  \code{alpha=1}, else it is the enet, mnet or snet sequence.
 For the other families, this is a lasso (mcp, scad) or elastic net (mnet, snet) regularization path
  for fitting the generalized linear regression
  paths, by maximizing the appropriate penalized log-likelihood.
  Note that the objective function for \code{"gaussian"} is \deqn{1/2*
  weights*RSS + \lambda*penalty,} if \code{standardize=FALSE} and \deqn{1/2*
  \frac{weights}{\sum(weights)}*RSS + \lambda*penalty,} if \code{standardize=TRUE}. For the other models it is
  \deqn{-\sum (weights * loglik) + \lambda*penalty} if \code{standardize=FALSE} and \deqn{-\frac{weights}{\sum(weights)} * loglik + \lambda*penalty} if \code{standardize=TRUE}. 
  }
\value{
An object with S3 class \code{"glmreg"} for the various types of models.
  \item{call}{the call that produced the model fit}
  \item{b0}{Intercept sequence of length \code{length(lambda)}}
  \item{beta}{A \code{nvars x length(lambda)} matrix of coefficients.}
  \item{lambda}{The actual sequence of \code{lambda} values used}
  \item{satu}{satu=1 if a saturated model (deviance/null deviance < 0.05) is fit. Otherwise satu=0. The number of \code{nlambda} sequence may be truncated before \code{nlambda} is reached if \code{satu=1}.}
  \item{dev}{The computed deviance (for \code{"gaussian"}, this
      is the R-square). The deviance calculations incorporate weights if
  present in the model. The deviance is defined to be 2*(loglike_sat -
  loglike), where loglike_sat is the log-likelihood for the saturated
  model (a model with a free parameter per observation). }
    \item{nulldev}{Null deviance (per observation). This is defined to
  be  2*(loglike_sat -loglike(Null)); The NULL model refers to the
  intercept model.}
  \item{nobs}{number of observations}
}
\references{
Breheny, P. and Huang, J. (2011) \emph{Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection}. \emph{Ann. Appl. Statist.}, \bold{5}: 232-253.

Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)
\emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]
}

\author{Zhu Wang <zwang@connecticutchildrens.org>}

\seealso{\code{\link{glmreg}}}
\keyword{models}
\keyword{regression}

 
