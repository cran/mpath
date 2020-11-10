\name{nclreg}
\alias{nclreg}
\alias{nclreg.formula}
\alias{nclreg.matrix}
\alias{nclreg.default}

\title{Optimize a nonconvex loss with regularization}
\description{
    Fit a linear model via penalized nonconvex loss function.  The
  regularization path is computed for the lasso (or elastic net penalty), scad (or snet) and mcp (or mnet penalty), at a grid
  of values for the regularization parameter lambda. The name refers to \bold{N}on\bold{C}onvex \bold{L}oss with \bold{REG}ularization.
  }
\usage{
\method{nclreg}{formula}(formula, data, weights, offset=NULL, contrasts=NULL, ...)
\method{nclreg}{matrix}(x, y, weights, offset=NULL, ...)
\method{nclreg}{default}(x,  ...)
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
  \item{...}{Other arguments passing to \code{nclreg_fit}}
}
\details{
  The sequence of robust models implied by \code{lambda} is fit by majorization-minimization along with coordinate
  descent. Note that the objective function is \deqn{weights*loss + \lambda*penalty,} if \code{standardize=FALSE} and \deqn{\frac{weights}{\sum(weights)}*loss + \lambda*penalty,} if \code{standardize=TRUE}. 
  }

\value{
An object with S3 class \code{"nclreg"} for the various types of models.
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
 Zhu Wang (2019)
 \emph{MM for Penalized Estimation},
\url{https://arxiv.org/abs/1912.11119 
}}

\author{Zhu Wang <wangz1@uthscsa.edu>}
\seealso{\code{\link{print}}, \code{\link{predict}}, \code{\link{coef}} and \code{\link{plot}} methods, and the \code{\link{cv.nclreg}} function.}
\examples{
#binomial
x=matrix(rnorm(100*20),100,20)
g2=sample(c(-1,1),100,replace=TRUE)
### different solution paths via a combination of type.path, decreasing and type.init
fit1=nclreg(x,g2,s=1,rfamily="closs",type.path="active",decreasing=TRUE,type.init="bst")
fit2=nclreg(x,g2,s=1,rfamily="closs",type.path="active",decreasing=FALSE,type.init="bst")
fit3=nclreg(x,g2,s=1,rfamily="closs",type.path="nonactive",decreasing=TRUE,type.init="bst")
fit4=nclreg(x,g2,s=1,rfamily="closs",type.path="nonactive",decreasing=FALSE,type.init="bst")
fit5=nclreg(x,g2,s=1,rfamily="closs",type.path="active",decreasing=TRUE,type.init="ncl")
fit6=nclreg(x,g2,s=1,rfamily="closs",type.path="active",decreasing=FALSE,type.init="ncl")
fit7=nclreg(x,g2,s=1,rfamily="closs",type.path="nonactive",decreasing=TRUE,type.init="ncl")
fit8=nclreg(x,g2,s=1,rfamily="closs",type.path="nonactive",decreasing=FALSE,type.init="ncl")
}
\keyword{models}
\keyword{regression}

 
