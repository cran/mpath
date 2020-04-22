\name{cv.zipath}
\alias{cv.zipath}
\alias{cv.zipath.matrix}
\alias{cv.zipath.formula}
\alias{cv.zipath.default}
\alias{predict.cv.zipath}
\alias{coef.cv.zipath}
\title{Cross-validation for zipath}
\description{Does k-fold cross-validation for zipath, produces a plot,
  and returns cross-validated log-likelihood values for \code{lambda}}
\usage{
\method{cv.zipath}{formula}(formula, data, weights, offset=NULL, contrasts=NULL, ...)
\method{cv.zipath}{matrix}(X, Z, Y, weights, offsetx=NULL, offsetz=NULL, ...)
\method{cv.zipath}{default}(X, ...)
\method{predict}{cv.zipath}(object, newdata, \dots)
\method{coef}{cv.zipath}(object, which=object$lambda.which, model = c("full", "count", "zero"), \dots)
}
\arguments{
  \item{formula}{symbolic description of the model with an optional numeric vector \code{offset} with an a priori known component to be included in the linear predictor of the count model or zero model. Offset must be a variable in \code{data} if used, while this is optional in \code{zipath}. See an example below.}
  \item{data}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
   \item{offset}{optional numeric vector with an a priori known component to be
     included in the linear predictor of the count model or zero model. See below   for an example.}
   \item{X}{predictor matrix of the count model}
    \item{Z}{predictor matrix of the zero model}
    \item{Y}{response variable}
   \item{offsetx, offsetz}{optional numeric vector with an a priori known component to be
     included in the linear predictor of the count model (offsetx)or zero model (offsetz).}
  \item{contrasts}{a list with elements \code{"count"} and \code{"zero"}
      containing the contrasts corresponding to \code{levels} from the
      respective models} 
 \item{object}{object of class \code{cv.zipath}.}
 \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
 \item{which}{Indices of the pair of penalty parameters \code{lambda.count} and \code{lambda.zero} at which
    estimates are extracted.  By default, the one which generates the optimal cross-validation value.}
 \item{model}{character specifying for which component of the model the
    estimated coefficients should be extracted.}
\item{\dots}{Other arguments that can be passed to \code{zipath}.}
}
\details{The function runs \code{zipath} \code{nfolds}+1 times; the
  first to compute the \code{(lambda.count, lambda.zero)} sequence, and then to
  compute the fit with each of the folds omitted. The log-likelihood value is
  accumulated, and the average value and standard deviation over the
  folds is computed. Note that \code{cv.zipath} can be used to search for
  values for \code{count.alpha} or \code{zero.alpha}: it is required to call \code{cv.zipath} with a fixed vector \code{foldid} for different values of \code{count.alpha} or \code{zero.alpha}.

The method for \code{\link[stats]{coef}} by default
  return a single vector of coefficients, i.e., all coefficients are concatenated. By setting the \code{model}
  argument, the estimates for the corresponding model components can be extracted.
}
\value{an object of class \code{"cv.zipath"} is returned, which is a
  list with the components of the cross-validation fit.
 \item{fit}{a fitted zipath object for the full data.}
 \item{residmat}{matrix for cross-validated log-likelihood at each \code{(count.lambda, zero.lambda)} sequence}
  \item{bic}{matrix of BIC values with row values for \code{lambda} and column values for \code{k}th cross-validation}
  \item{cv}{The mean cross-validated log-likelihood - a vector of length
    \code{length(count.lambda)}.}
  \item{cv.error}{estimate of standard error of \code{cv}.}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in.} 
 \item{lambda.which}{index of \code{(count.lambda, zero.lambda)} that gives maximum \code{cv}.}
 \item{lambda.optim}{value of \code{(count.lambda, zero.lambda)} that gives maximum \code{cv}.}
}
\references{
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)      \emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]

Zhu Wang, Shuangge Ma, Ching-Yun Wang, Michael Zappitelli, Prasad Devarajan and Chirag R. Parikh (2014)
\emph{EM for Regularized Zero Inflated Regression Models with Applications to Postoperative Morbidity after Cardiac Surgery in Children}, \emph{Statistics in Medicine}. 33(29):5192-208.

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}

\author{Zhu Wang <wangz1@uthscsa.edu>}

\seealso{\code{\link{zipath}} and \code{\link{plot}}, \code{\link{predict}}, and \code{\link{coef}} methods for \code{"cv.zipath"} object.}
\examples{
\dontrun{
data("bioChemists", package = "pscl")
fm_zip <- zipath(art ~ . | ., data = bioChemists, family = "poisson", nlambda=10)
fm_cvzip <- cv.zipath(art ~ . | ., data = bioChemists, family = "poisson", nlambda=10)
### prediction from the best model
pred <- predict(fm_zip, newdata=bioChemists, which=fm_cvzip$lambda.which)
coef(fm_zip, which=fm_cvzip$lambda.which)
fm_znb <- zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10)
fm_cvznb <- cv.zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10)
pred <- predict(fm_znb, which=fm_cvznb$lambda.which)
coef(fm_znb, which=fm_cvznb$lambda.which)
fm_zinb2 <- zipath(art ~ . +offset(log(phd))| ., data = bioChemists, 
		      family = "poisson", nlambda=10)
fm_cvzinb2 <- cv.zipath(art ~ . +offset(log(phd))| ., data = bioChemists, 
		      family = "poisson", nlambda=10)
pred <- predict(fm_zinb2, which=fm_cvzinb2$lambda.which)
coef(fm_zinb2, which=fm_cvzinb2$lambda.which)
}
}
\keyword{models}
\keyword{regression}

