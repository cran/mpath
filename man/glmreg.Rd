\name{glmreg}
\alias{glmreg}
\alias{glmreg.formula}
\alias{glmreg.matrix}
\alias{glmreg.default}
\alias{logLik.glmreg}
\alias{deviance.glmreg}

\title{fit a GLM with lasso (or elastic net), snet or mnet regularization}
\description{
    Fit a generalized linear model via penalized maximum likelihood.  The
  regularization path is computed for the lasso (or elastic net penalty), scad (or snet) and mcp (or mnet penalty), at a grid
  of values for the regularization parameter lambda. Fits
  linear, logistic, Poisson and negative binomial (fixed scale parameter) regression models.
  }
\usage{
\method{glmreg}{formula}(formula, data, weights, offset=NULL, contrasts=NULL, 
x.keep=FALSE, y.keep=TRUE, ...)
\method{glmreg}{matrix}(x, y, weights, offset=NULL, ...)
\method{glmreg}{default}(x,  ...)
}

\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights. If \code{standardize=TRUE}, weights are renormalized to weights/sum(weights). If \code{standardize=FALSE}, weights are kept as original input}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of cases. Currently only one offset term can be included in the formula.}
  \item{x}{input matrix, of dimension nobs x nvars; each row is an
  observation vector}
  \item{y}{response variable. Quantitative for \code{family="gaussian"}.
  Non-negative counts for \code{family="poisson"} or \code{family="negbin"}. For
  \code{family="binomial"} should be either a factor with two levels or
  a vector of proportions.}
  \item{x.keep, y.keep}{logical values: keep response variables or keep response variable?} 
  \item{contrasts}{the contrasts corresponding to \code{levels} from the
    respective models}
  \item{...}{Other arguments passing to \code{glmreg_fit}}
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
  \item{call}{the call that produced this object}
  \item{b0}{Intercept sequence of length \code{length(lambda)}}
  \item{beta}{A \code{nvars x
      length(lambda)} matrix of coefficients.}
  \item{lambda}{The actual sequence of \code{lambda} values used}
  \item{dev}{The computed deviance (for \code{"gaussian"}, this
      is the R-square). The deviance calculations incorporate weights if
  present in the model. The deviance is defined to be 2*(loglike_sat -
  loglike), where loglike_sat is the log-likelihood for the saturated
  model (a model with a free parameter per observation). }
    \item{nulldev}{Null deviance (per observation). This is defined to
  be  2*(loglike_sat -loglike(Null)); The NULL model refers to the
  intercept model.}
  \item{nobs}{number of observations}
  \item{pll}{penalized log-likelihood values for standardized coefficients in the IRLS iterations. For \code{family="gaussian"}, not implemented yet.}
  \item{pllres}{penalized log-likelihood value for the estimated model on the original scale of coefficients}
  \item{fitted.values}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function. }
}
\references{
Breheny, P. and Huang, J. (2011) \emph{Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection}. \emph{Ann. Appl. Statist.}, \bold{5}: 232-253.

Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014) 
\emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}.
2014 Apr 17. [Epub ahead of print]
}

\author{Zhu Wang <wangz1@uthscsa.edu>}
\seealso{\code{\link{print}}, \code{\link{predict}}, \code{\link{coef}} and \code{\link{plot}} methods, and the \code{\link{cv.glmreg}} function.}
\examples{
#binomial
x=matrix(rnorm(100*20),100,20)
g2=sample(0:1,100,replace=TRUE)
fit2=glmreg(x,g2,family="binomial")
#poisson and negative binomial
data("bioChemists", package = "pscl")
fm_pois <- glmreg(art ~ ., data = bioChemists, family = "poisson")
coef(fm_pois)
fm_nb1 <- glmreg(art ~ ., data = bioChemists, family = "negbin", theta=1)
coef(fm_nb1)
#offset
x <- matrix(rnorm(100*20),100,20)
y <- rpois(100, lambda=1)
exposure <- rep(0.5, length(y))
fit2 <- glmreg(x,y, lambda=NULL, nlambda=10, lambda.min.ratio=1e-4, 
	       offset=log(exposure), family="poisson")
predict(fit2, newx=x, newoffset=log(exposure))
\dontrun{
fm_nb2 <- glmregNB(art ~ ., data = bioChemists)
coef(fm_nb2)
}
}
\keyword{models}
\keyword{regression}

 
