\name{ccglm}
\alias{ccglm}
\alias{ccglm.formula}

\title{fit a CC-estimator for robust generalized linear models}
\description{
    Fit a CC-estimator for robust generalized linear models
  }
\usage{
\method{ccglm}{formula}(formula, data, weights, offset=NULL, contrasts=NULL,
 cfun="ccave", dfun=gaussian(), s=NULL, delta=0.1, fk=NULL, init.family=NULL,
 iter=10, reltol=1e-5, theta, x.keep=FALSE, y.keep=TRUE, trace=FALSE, ...)
}

\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights.}
  \item{x}{input matrix, of dimension nobs x nvars; each row is an
  observation vector}
  \item{y}{response variable. Quantitative for \code{dfun=1} and -1/1 for classification. }
  \item{contrasts}{the contrasts corresponding to \code{levels} from the
    respective models}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of cases. Currently only one offset term can be included in the formula.}
   \item{cfun}{ character, type of convex cap (concave) function.\cr
   Valid options are:
   \itemize{
       \item \code{"hcave"}
       \item \code{"acave"}
       \item \code{"bcave"}
       \item \code{"ccave"}
       \item \code{"dcave"}
       \item \code{"ecave"}
       \item \code{"gcave"}
       \item \code{"tcave"}
   }
   }
   \item{dfun}{ character, type of convex component. \cr
   Valid options are:
   \itemize{
       \item \code{gaussian()}
       \item \code{binomial()}
       \item \code{poisson()}
       }
   }
   \item{init.family}{ character value for initial family, one of "clossR","closs","gloss","qloss", which can be used to derive an initial estimator, if the selection is different from the default value }
  \item{s}{ tuning parameter of \code{cfun}. \code{s > 0} and can be equal to 0 for \code{cfun="tcave"}. If \code{s} is too close to 0 for  \code{cfun="acave", "bcave", "ccave"}, the calculated weights can become 0 for all observations, thus crash the program.} 
 \item{delta}{a small positive number provided by user only if \code{cfun="gcave"} and \code{0 < s <1}}
  \item{fk}{ predicted values at an iteration in the COCO algorithm }
 \item{iter}{number of iteration in the COCO algorithm}
 \item{reltol}{convergency criteria in the COCO algorithm}
\item{theta}{ an overdispersion scaling parameter for \code{family=negbin()}} 
  \item{x.keep, y.keep}{ logical values indicating whether the response
          vector and model matrix used in the fitting process should be
          returned as components of the returned value, x is a design matrix of dimension n * p,
          and x is a vector of observations of length n.}
 \item{trace}{if \code{TRUE}, fitting progress is reported}
  \item{...}{other arguments passing to \code{ccglm}}
}
\details{
  A robust linear, logistic or Poisson regression model is fit by the COCO algorithm, where the loss function is a composite function of \code{cfun}o\code{dfun}.
  }
\value{
An object with S3 class \code{"ccglm", "glm"} for various types of models.
  \item{call}{the call that produced the model fit}
  \item{weights}{original weights used in the model}
  \item{weights_update}{weights in the final iteration of the COCO algorithm}
  \item{cfun, s, dfun}{original input arguments}
  \item{is.offset}{is offset used?}
}

\references{
  Zhu Wang (2020)
  \emph{Unified Robust Estimation via the COCO, arXiv e-prints}, \url{https://arxiv.org/abs/2010.02848}
 }

\author{Zhu Wang <wangz1@uthscsa.edu>}
\seealso{\code{\link{print}}, \code{\link{predict}}, \code{\link{coef}}.}
\examples{
x=matrix(rnorm(100*20),100,20)
g2=sample(c(-1,1),100,replace=TRUE)
fit=ccglm(g2~x,data=data.frame(cbind(x, g2)), s=1, cfun="ccave", dfun=gaussian())
fit$weights_update
}
 
\keyword{models}
\keyword{regression}

 
