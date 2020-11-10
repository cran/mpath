\name{ccsvm}
\alias{ccsvm}
\alias{ccsvm.formula}
\alias{ccsvm.matrix}
\alias{ccsvm.default}
\alias{coef.ccsvm}

\title{fit case weighted support vector machines with robust loss functions}
\description{
    Fit case weighted support vector machines with robust loss functions.
  }
\usage{
\method{ccsvm}{formula}(formula, data, weights, contrasts=NULL, ...)
\method{ccsvm}{matrix}(x, y, weights, ...)
\method{ccsvm}{default}(x,  ...)
}

\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights}
  \item{x}{input matrix, of dimension nobs x nvars; each row is an
  observation vector}
  \item{y}{response variable. Quantitative for \code{type="eps-regression", "nu-regression"} and -1/1 for \code{type="C-classification", "nu-Classification"}.}
   \item{contrasts}{the contrasts corresponding to \code{levels} from the
     respective models}
  \item{...}{Other arguments passing to \code{ccsvm_fit}}
}
\details{
  The model is fit by the COCO algorithm.

  For linear kernel, the coefficients of the regression/decision hyperplane
   can be extracted using the \code{coef} method.
}

\value{
An object with S3 class \code{"wsvm"} for various types of models.
\item{call}{the call that produced this object}
\item{weights_update}{weights in the final iteration of the COCO algorithm}
\item{cfun, s}{original input arguments}
\item{delta}{delta value used for \code{cfun="gcave"}}
}

\references{
   Zhu Wang (2020)
   \emph{Unified Robust Estimation via the COCO, arXiv e-prints}, \url{https://arxiv.org/abs/2010.02848}
  }

\author{Zhu Wang <wangz1@uthscsa.edu>}
\seealso{\code{\link{print}}, \code{\link{predict}}, \code{\link{coef}}.}
\examples{
#binomial
x=matrix(rnorm(100*20),100,20)
g2=sample(c(-1,1),100,replace=TRUE)
fit=ccsvm(x,g2,s=1,cfun="ccave",type="C-classification")
}
\keyword{models}
\keyword{regression}

 
