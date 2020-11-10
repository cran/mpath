\name{loss3}
\alias{loss3}
\title{Composite Loss Value for GLM}
\description{
 Compute composite loss value
}
\usage{
loss3(y, mu, theta, weights, cfun, family, s, delta)
}
\arguments{
 \item{y}{response variable values, 0/1 if \code{family=2}, or binomial}
 \item{mu}{response prediction of \code{y}. If \code{mu} is linear predictor, use function \code{loss2} instead}
 \item{theta}{scale parameter for \code{family=4}, negative binomial}
 \item{weights}{observation weights, same length as \code{y}}
 \item{cfun}{integer from 1-8, concave function as in \code{ccglm_fit}}
 \item{family}{integer 2, 3 or 4, convex function binomial, Poisson or negative binomial, respectively}
\item{s}{ tuning parameter of \code{cfun}. \code{s > 0} and can be equal to 0 for \code{cfun="tcave"}.} 
   \item{delta}{a small positive number provided by user only if \code{cfun="gcave"} and \code{0 < s <1}}
}

\value{Weighted loss values}
\details{For large \code{s} values, the loss can be 0 with \code{cfun=2,3,4}, or "acave", "bcave", "ccave".}
\seealso{
    \code{\link{loss2}}
    \code{\link{ccglm}}
    \code{\link{ccglmreg}}
    \code{\link{loss2_ccsvm}}
}
\references{
   Zhu Wang (2020)
   \emph{Unified Robust Estimation via the COCO, arXiv e-prints}, \url{https://arxiv.org/abs/2010.02848}
  }
\author{Zhu Wang <wangz1@uthscsa.edu>}

\keyword{regression}
\keyword{classification}
