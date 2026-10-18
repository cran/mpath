\name{irsvm}
\alias{irsvm}
\alias{irsvm.formula}
\alias{irsvm.matrix}
\alias{irsvm.default}
\alias{coef.irsvm}

\title{fit case weighted support vector machines with robust loss functions}
\description{
    Fit case weighted support vector machines with robust loss functions. This is the wrapper function of \code{irsvm_fit}, which does the computing.
  }
\usage{
\method{irsvm}{formula}(formula, data, weights, contrasts=NULL, ...)
\method{irsvm}{matrix}(x, y, weights, ...)
\method{irsvm}{default}(x,  ...)
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
  \item{...}{Other arguments passing to \code{irsvm_fit}}
}
\details{
  Fit a robust SVM where the loss function is a composite function \code{cfun}o\code{type} + penalty.
    The model is fit by the iteratively reweighted SVM, an application of the iteratively reweighted convex optimization (IRCO). Here convex is the loss function  induced by \code{type}. 

  For linear kernel, the coefficients of the regression/decision hyperplane
   can be extracted using the \code{coef} method.
}

\value{
An object with S3 class \code{"wsvm"} for various types of models.
\item{call}{the call that produced this object}
\item{weights_update}{weights in the final iteration of the IRCO algorithm}
\item{cfun, s}{original input arguments}
\item{delta}{delta value used for \code{cfun="gcave"}}
}

\references{
   Zhu Wang (2024)
    \emph{Unified Robust Estimation}, \emph{Australian & New Zealand Journal of      Statistics}. 66(1):77-102.
}

\author{Zhu Wang <zwang145@uthsc.edu>}
\seealso{\code{\link{irsvm_fit}}, \code{\link{print}}, \code{\link{predict}}, \code{\link{coef}}.}
\examples{
#binomial
x=matrix(rnorm(100*20),100,20)
g2=sample(c(-1,1),100,replace=TRUE)
fit=irsvm(x,g2,s=1,cfun="ccave",type="C-classification")
}
\keyword{models}
\keyword{regression}

 
