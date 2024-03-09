\name{irglmreg}
\alias{irglmreg}
\alias{irglmreg.formula}
\alias{irglmreg.matrix}
\alias{irglmreg.default}

\title{Fit a robust penalized generalized linear models}
\description{
    Fit a robust penalized GLM where the loss function is a composite function \code{cfun}o\code{dfun} + penalty. This is the wrapper function of \code{irglmreg_fit} 
  }
\usage{
\method{irglmreg}{formula}(formula, data, weights, offset=NULL, contrasts=NULL, ...)
\method{irglmreg}{matrix}(x, y, weights, offset=NULL, ...)
\method{irglmreg}{default}(x,  ...)
}

\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights. If \code{standardize=TRUE}, weights are renormalized to weights/sum(weights). If \code{standardize=FALSE}, weights are kept as original input}
  \item{x}{input matrix, of dimension nobs x nvars; each row is an
  observation vector}
  \item{y}{response variable. Quantitative for \code{rfamily="clossR"} and -1/1 for classification. }
  \item{offset}{Not implemented yet}
  \item{contrasts}{the contrasts corresponding to \code{levels} from the
    respective models}
  \item{...}{Other arguments passing to \code{irglmreg_fit}}
}
\details{
    The computing is done by the iteratively reweighted penalized GLM, an application of the iteratively reweighted convex optimization (IRCO). Here convex is the loss function induced by \code{dfun}, not the penalty function. The output \code{weights_update} is a useful diagnostic to the outlier status of     the observations.
The regularization path is computed for the lasso (or elastic net penalty), scad (or snet) and mcp (or mnet penalty), at a grid of values for the regularization parameter lambda. 
  The sequence of robust models implied by \code{lambda} is fit by the IRCO along with coordinate
  descent. Note that the objective function is \deqn{weights*loss + \lambda*penalty,} if \code{standardize=FALSE} and \deqn{
  \frac{weights}{\sum(weights)}*loss + \lambda*penalty,} if \code{standardize=TRUE}. 
  }

\value{
An object with S3 class \code{"irglmreg"} for the various types of models.
  \item{call}{the call that produced this object}
  \item{b0}{Intercept sequence of length \code{length(lambda)}}
  \item{beta}{A \code{nvars x
      length(lambda)} matrix of coefficients.}
  \item{lambda}{The actual sequence of \code{lambda} values used}
  \item{nobs}{number of observations}
  \item{risk}{if \code{type.path="nonactive"}, a matrix with number of rows \code{iter} and number of columns \code{nlambda}, loss values along the regularization path. If \code{type.path="fast"}, a vector of length \code{nlambda}, loss values along the regularization path}
  \item{pll}{if \code{type.path="nonactive"}, a matrix with number of rows \code{iter} and number of columns \code{nlambda}, penalized loss values along the regularization path. If \code{type.path="fast"}, a vector of length \code{nlambda}, penalized loss values along the regularization path}
  \item{fitted.values}{predicted values depending on \code{standardize}, internal use only}
}

\references{
   Zhu Wang (2024)
    \emph{Unified Robust Estimation}, \emph{Australian & New Zealand Journal of      Statistics}. 66(1):77-102.
}

\author{Zhu Wang <zwang145@uthsc.edu>}
\seealso{\code{\link{print}}, \code{\link{predict}}, \code{\link{coef}} and \code{\link{plot}} methods, and the \code{\link{cv.irglmreg}} function.}
\examples{
#binomial
x=matrix(rnorm(100*20),100,20)
g2=sample(c(-1,1),100,replace=TRUE)
fit1=irglmreg(x,g2,s=1,cfun="ccave",dfun="gaussian",type.path="active",
              decreasing=TRUE,type.init="bst")
#fit1$risk
\dontrun{
### different solution paths via a combination of type.path, decreasing and type.init
fit1=irglmreg(x,g2,s=1,cfun="ccave",dfun="gaussian",type.path="active",
           decreasing=TRUE,type.init="bst")
fit2=irglmreg(x,g2,s=1,cfun="ccave",dfun="gaussian",type.path="active",
           decreasing=FALSE,type.init="bst")
fit3=irglmreg(x,g2,s=1,cfun="ccave",dfun="gaussian",type.path="nonactive",
           decreasing=TRUE,type.init="bst")
fit4=irglmreg(x,g2,s=1,cfun="ccave",dfun="gaussian",type.path="nonactive",
           decreasing=FALSE,type.init="bst")
fit5=irglmreg(x,g2,s=1,cfun="ccave",dfun="gaussian",type.path="active",
           decreasing=TRUE,type.init="co")
fit6=irglmreg(x,g2,s=1,cfun="ccave",dfun="gaussian",type.path="active",
           decreasing=FALSE,type.init="co")
fit7=irglmreg(x,g2,s=1,cfun="ccave",dfun="gaussian",type.path="nonactive",
           decreasing=TRUE,type.init="co")
fit8=irglmreg(x,g2,s=1,cfun="ccave",dfun="gaussian",type.path="nonactive",
           decreasing=FALSE,type.init="co")
}
}
\keyword{models}
\keyword{regression}

 
