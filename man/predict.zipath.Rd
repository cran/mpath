\name{predict.zipath}
\alias{predict.zipath}
\alias{residuals.zipath}
\alias{terms.zipath}
\alias{model.matrix.zipath}
\alias{coef.zipath}
\alias{summary.zipath}
\alias{print.summary.zipath}
\alias{logLik.zipath}
\alias{fitted.zipath}
\alias{predprob.zipath}

\title{Methods for zipath Objects}
\description{
  Methods for extracting information from fitted penalized zero-inflated
  regression model objects of class \code{"zipath"}.
}

\usage{
\method{predict}{zipath}(object, newdata,  which = 1:object$nlambda, 
  type = c("response", "prob", "count", "zero", "nonzero"), na.action = na.pass,
  at = NULL, \dots)
\method{residuals}{zipath}(object, type = c("pearson", "response"), \dots)

\method{coef}{zipath}(object, which=1:object$nlambda, model = c("full", "count", "zero"), \dots)

\method{terms}{zipath}(x, model = c("count", "zero"), \dots)
\method{model.matrix}{zipath}(object, model = c("count", "zero"), \dots)
}

\arguments{
  \item{object, x}{an object of class \code{"zipath"} as returned by
    \code{\link{zipath}}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{which}{Indices of the penalty parameters \code{lambda} at which
    predictions are required.  By default, all indices are returned.}
  \item{type}{character specifying the type of predictions or residuals,
    respectively. For details see below.}
  \item{na.action}{function determining what should be done with missing values
    in \code{newdata}. The default is to predict \code{NA}.}
  \item{at}{optionally, if \code{type = "prob"}, a numeric vector at which
    the probabilities are evaluated. By default \code{0:max(y)} is used
    where \code{y} is the original observed response.}
  \item{model}{character specifying for which component of the model the
    terms or model matrix should be extracted.}
  \item{\dots}{currently not used.}
}

\details{
  Re-uses the design of function zeroinfl in package pscl (see reference). A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"zipath"}, including methods to the generic functions
  \code{\link[base]{print}} and \code{\link[base]{summary}} which print the estimated
  coefficients along with some further information. 
  As usual, the \code{summary} method returns an object of class \code{"summary.zipath"}
  containing the relevant summary statistics which can subsequently be printed
  using the associated \code{print} method.
  
  The methods for \code{\link[stats]{coef}} by default
  return a single vector of coefficients and their associated covariance matrix,
  respectively, i.e., all coefficients are concatenated. By setting the \code{model}
  argument, the estimates for the corresponding model components can be extracted.
  
  Both the \code{\link[stats]{fitted}} and \code{\link[stats]{predict}} methods can
  compute fitted responses. The latter additionally provides the predicted density
  (i.e., probabilities for the observed counts), the predicted mean from the count
  component (without zero inflation) and the predicted probability for the zero
  component. The \code{\link[stats]{residuals}} method can compute
  raw residuals (observed - fitted) and Pearson residuals (raw residuals scaled by
  square root of variance function).
}

\references{
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)      \emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]

Zhu Wang, Shuangge Ma, Ching-Yun Wang, Michael Zappitelli, Prasad Devarajan and Chirag R. Parikh (2014)
\emph{EM for Regularized Zero Inflated Regression Models with Applications to Postoperative Morbidity after Cardiac Surgery in Children}, \emph{Statistics in Medicine}. 33(29):5192-208.

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}
\author{Zhu Wang <wangz1@uthscsa.edu>}

\seealso{\code{\link{zipath}}}

\examples{
\dontrun{
data("bioChemists", package = "pscl")
fm_zip <- zipath(art ~ . | ., data = bioChemists, nlambda=10)
plot(residuals(fm_zip) ~ fitted(fm_zip))
coef(fm_zip, model = "count")
coef(fm_zip, model = "zero")
summary(fm_zip)
logLik(fm_zip)
}
}

\keyword{regression}
