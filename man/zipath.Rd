\name{zipath}
\alias{zipath}
\alias{zipath.formula}
\alias{zipath.matrix}
\alias{zipath.default}

\title{Fit zero-inflated count data linear model with lasso (or elastic net), snet or mnet regularization}
\description{
  Fit zero-inflated regression models for count data via penalized maximum likelihood.
}

\usage{
\method{zipath}{formula}(formula, data, weights, offset=NULL, contrasts=NULL, ... )
\method{zipath}{matrix}(X, Z, Y, weights, offsetx=NULL, offsetz=NULL, ...)
\method{zipath}{default}(X, ...)
}
\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights.}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor of the count model or zero model. See below for an example.}
  \item{contrasts}{a list with elements \code{"count"} and \code{"zero"}
     containing the contrasts corresponding to \code{levels} from the
     respective models}
\item{X}{predictor matrix of the count model}
   \item{Z}{predictor matrix of the zero model}
   \item{Y}{response variable}
  \item{offsetx, offsetz}{optional numeric vector with an a priori known component to be
    included in the linear predictor of the count model (offsetx)or zero model (offsetz).}
\item{...}{Other arguments which can be passed to \code{glmreg} or \code{glmregNB}}
}

\value{
  An object of class \code{"zipath"}, i.e., a list with components including
  \item{coefficients}{a list with elements \code{"count"} and \code{"zero"}
    containing the coefficients from the respective models,}
  \item{residuals}{a vector of raw residuals (observed - fitted),}
  \item{fitted.values}{a vector of fitted means,}
  \item{weights}{the case weights used,}
  \item{terms}{a list with elements \code{"count"}, \code{"zero"} and
    \code{"full"} containing the terms objects for the respective models,}
  \item{theta}{estimate of the additional \eqn{\theta}{theta} parameter of the
    negative binomial model (if a negative binomial regression is used),}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{family}{character string describing the count distribution used,}
  \item{link}{character string describing the link of the zero-inflation model,}
  \item{linkinv}{the inverse link function corresponding to \code{link},}
  \item{converged}{logical value, TRUE indicating successful convergence of \code{zipath}, FALSE indicating otherwise}
  \item{call}{the original function call}
  \item{formula}{the original formula}
  \item{levels}{levels of the categorical regressors}
  \item{contrasts}{a list with elements \code{"count"} and \code{"zero"}
    containing the contrasts corresponding to \code{levels} from the
    respective models,}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the response count vector (if \code{y = TRUE}),}
  \item{x}{a list with elements \code{"count"} and \code{"zero"}
    containing the model matrices from the respective models
    (if \code{x = TRUE}),}
}

\references{     
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)      \emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]

Zhu Wang, Shuangge Ma, Ching-Yun Wang, Michael Zappitelli, Prasad Devarajan and Chirag R. Parikh (2014)
\emph{EM for Regularized Zero Inflated Regression Models with Applications to Postoperative Morbidity after Cardiac Surgery in Children}, \emph{Statistics in Medicine}. 33(29):5192-208.

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}

\author{Zhu Wang <wangz1@uthscsa.edu>}

\seealso{\code{\link{zipath_fit}},
  \code{\link{glmreg}}, \code{\link{glmregNB}}
}

\examples{
## data
data("bioChemists", package = "pscl")
## with simple inflation (no regressors for zero component)
fm_zip <- zipath(art ~ 1 | ., data = bioChemists, nlambda=10)
summary(fm_zip)
fm_zip <- zipath(art ~ . | 1, data = bioChemists, nlambda=10)
summary(fm_zip)
\dontrun{
fm_zip <- zipath(art ~ . | 1, data = bioChemists, nlambda=10)
summary(fm_zip)
fm_zinb <- zipath(art ~ . | 1, data = bioChemists, family = "negbin", nlambda=10)
summary(fm_zinb)
## inflation with regressors
## ("art ~ . | ." is "art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment")
fm_zip2 <- zipath(art ~ . | ., data = bioChemists, nlambda=10)
summary(fm_zip2)
fm_zinb2 <- zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10)
summary(fm_zinb2)
### non-penalized regression, compare with zeroinfl
fm_zinb3 <- zipath(art ~ . | ., data = bioChemists, family = "negbin", 
lambda.count=0, lambda.zero=0, reltol=1e-12)
summary(fm_zinb3)
library("pscl")
fm_zinb4 <- zeroinfl(art ~ . | ., data = bioChemists, dist = "negbin")
summary(fm_zinb4)
### offset
exposure <- rep(0.5, dim(bioChemists)[1])
fm_zinb <- zipath(art ~ . +offset(log(exposure))| ., data = bioChemists, 
		  family = "poisson", nlambda=10)
coef <- coef(fm_zinb)
### offset can't be specified in predict function as it has been contained
pred <- predict(fm_zinb)
## without inflation
## ("art ~ ." is "art ~ fem + mar + kid5 + phd + ment")
fm_pois <- glmreg(art ~ ., data = bioChemists, family = "poisson")
coef <- coef(fm_pois)
fm_nb <- glmregNB(art ~ ., data = bioChemists)
coef <- coef(fm_nb)
### high-dimensional
bioChemists <- cbind(matrix(rnorm(915*100), nrow=915), bioChemists)
fm_zinb <- zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10)
}
}

\keyword{regression}
