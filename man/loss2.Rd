\name{loss2}
\alias{loss2}
\title{Composite Loss Value}
\description{
 Compute composite loss value
}
\usage{
loss2(y, f, weights, cfun, dfun, s, delta=0.0001)
}
\arguments{
 \item{y}{response variable values}
 \item{f}{linear predictor values of \code{y}. If \code{f} is predicted \code{response} of model, use function \code{loss3} instead}
 \item{weights}{observation weights, same length as \code{y}}
 \item{cfun}{integer from 1-8, concave function as in \code{ccglm_fit}}
 \item{dfun}{integer from 1-7, convex function as in \code{ccglm_fit}}
\item{s}{ tuning parameter of \code{cfun}. \code{s > 0} and can be equal    to 0 for \code{cfun="tcave"}.} 
   \item{delta}{a small positive number provided by user only if \code{cfun="gcave"} and \code{0 < s <1}}
}

\value{Weighted loss values}
\details{An internal function. For large \code{s} values, the loss can be 0 with \code{cfun=2,3,4}, or "acave", "bcave", "ccave".}
\seealso{
    \code{\link{loss3}}
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
