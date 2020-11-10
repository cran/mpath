\name{compute_wt}
\alias{compute_wt}
\title{Weight value from concave function}
\description{
 Weight value from concave function
}
\usage{
compute_wt(z, weights, cfun, s, delta=0.0001)
}
\arguments{
 \item{z}{vector nonnegative values from \code{dfun}, e.g., u^2/2}
 \item{weights}{optional numeric vector of weights.}
 \item{cfun}{integer from 1-8, concave function as in \code{ccglm_fit}}
 \item{s}{a numeric value, see details in \code{ccglm_fit}}
 \item{delta}{a positive small value, see details in \code{ccglm_fit}}
}

\value{Weight value from concave function}

\seealso{\code{\link{ccglmreg}}}
\examples{
u <- seq(0, 4, by=0.01)
z <- u^2/2 ### this is dfun
res <- compute_wt(z, cfun=1, s=1)
plot(z, res, ylab="Weight", type="l", lwd=2, 
     main=expression(paste("hcave", "   (", sigma, "=1)", )))
}
\references{
   Zhu Wang (2020)
   \emph{Unified Robust Estimation via the COCO, arXiv e-prints}, \url{https://arxiv.org/abs/2010.02848}
  }
\author{Zhu Wang <wangz1@uthscsa.edu>}

\keyword{regression}
\keyword{classification}
