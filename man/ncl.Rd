\name{ncl}
\alias{ncl}
\alias{ncl.formula}
\alias{ncl.matrix}
\alias{ncl.default}

\title{fit a nonconvex loss based robust linear model}
\description{
    Fit a linear model via penalized nonconvex loss function.  
  }
\usage{
\method{ncl}{formula}(formula, data, weights, offset=NULL, contrasts=NULL, 
x.keep=FALSE, y.keep=TRUE, ...)
\method{ncl}{matrix}(x, y, weights, offset=NULL, ...)
\method{ncl}{default}(x,  ...)
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
  \item{x.keep, y.keep}{ For glmreg: logical values indicating whether the response
          vector and model matrix used in the fitting process should be
          returned as components of the returned value.
          For ncl_fit: x is a design matrix of dimension n * p,
          and x is a vector of observations of length n.}
  \item{...}{Other arguments passing to \code{ncl_fit}}
}
\details{
  The robust linear model is fit by majorization-minimization along with linear regression. Note that the objective function is \deqn{weights*loss}.
  }

\value{
An object with S3 class \code{"ncl"} for the various types of models.
  \item{call}{the call that produced this object}
  \item{fitted.values}{predicted values}
  \item{h}{ pseudo response values in the MM algorithm }
}

\references{
   Zhu Wang (2021),
   \emph{MM for Penalized Estimation},
   \emph{TEST},
   \doi{10.1007/s11749-021-00770-2}
}

\author{Zhu Wang <wangz1@uthscsa.edu>}
\seealso{\code{\link{print}}, \code{\link{predict}}, \code{\link{coef}}.}
\examples{
#binomial
x=matrix(rnorm(100*20),100,20)
g2=sample(c(-1,1),100,replace=TRUE)
fit=ncl(x,g2,s=1,rfamily="closs")
}
\keyword{models}
\keyword{regression}

 
