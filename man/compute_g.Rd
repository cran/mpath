\name{compute_g}
\alias{compute_g}
\title{Compute concave function values}
\description{
 Compute concave function values
}
\usage{
compute_g(z, cfun, s, delta=0.0001)
}
\arguments{
 \item{z}{vector nonnegative values from \code{dfun}, e.g., u^2/2}
 \item{cfun}{integer from 1-8, concave function as in \code{irglm_fit}}
 \item{s}{a numeric value, see details in \code{irglmreg_fit}}
 \item{delta}{a positive small value, see details in \code{irglmreg_fit}}
}

\value{Concave function values}

\seealso{\code{\link{irglmreg}}}
\examples{
u <- seq(0, 4, by=0.01)
z <- u^2/2 ### this is dfun
res <- compute_g(z, cfun=1, s=1)
plot(z, res, ylab="Weight", type="l", lwd=2, 
     main=expression(paste("hcave", "   (", sigma, "=1)", )))
}
\references{
Zhu Wang (2024)
      \emph{Unified Robust Estimation}, \emph{Australian & New Zealand           Journal         of      Statistics}. 66(1):77-102.
}
\author{Zhu Wang <zwang145@uthsc.edu>}

\keyword{regression}
\keyword{classification}
