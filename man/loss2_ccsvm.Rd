\name{loss2_ccsvm}
\alias{loss2_ccsvm}
\title{Composite Loss Value for epsilon-insensitive Type}
\description{
 Compute composite loss value for epsilon-insensitive type function
}
\usage{
loss2_ccsvm(y, f, weights, cfun, dfun, s, eps, delta=0.0001)
}
\arguments{
 \item{y}{response variable values}
 \item{f}{fitted values of \code{y}}
 \item{weights}{observation weights, same length as \code{y}}
 \item{cfun}{integer from 1-8, concave function as in \code{ccsvm_fit}}
 \item{dfun}{integer value, only \code{dfun=2} is implemented for now. Convex function as in \code{ccsvm_fit}}
  \item{s}{ tuning parameter of \code{cfun}. \code{s > 0} and can be equal to 0 for \code{cfun="tcave"}.}
    \item{delta}{a small positive number provided by user only if  \code{cfun="gcave"} and \code{0 < s <1}}
 \item{eps}{non-negative parameter for epsilon-insensitive loss}
}

\value{Weighted loss values}
\details{For large \code{s} values, the loss can be 0 with \code{cfun=2,3,4},   or "acave", "bcave", "ccave".}
\seealso{
    \code{\link{ccglmreg}},
    \code{\link{loss2}}
}
\references{
   Zhu Wang (2020)
   \emph{Unified Robust Estimation via the COCO, arXiv e-prints}, \url{https://arxiv.org/abs/2010.02848}
  }
\author{Zhu Wang <wangz1@uthscsa.edu>}

\keyword{regression}
\keyword{classification}
