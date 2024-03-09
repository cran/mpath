\name{update_wt}
\alias{update_wt}
\title{Compute weight value}
\description{
 Compute weight value
}
\usage{
update_wt(y, ypre, weights, cfun, s, dfun, delta=0.0001)
}
\arguments{
 \item{y}{input value of response variable}
 \item{ypre}{predicted value of response variable}
 \item{weights}{optional numeric vector of weights.}
 \item{cfun}{integer from 1-8, concave function as in \code{irglm_fit}}
 \item{dfun}{integer value, convex function as in \code{irglm_fit}}
 \item{s}{a numeric value, see details in \code{irglm_fit}}
 \item{delta}{a positive small value, see details in \code{irglm_fit}}
}

\value{Weight value}

\seealso{\code{\link{compute_wt}}}
\references{
     Zhu Wang (2024)
    \emph{Unified Robust Estimation}, \emph{Australian & New Zealand Journal of      Statistics}. 66(1):77-102.
}
\author{Zhu Wang <zwang145@uthsc.edu>}

\keyword{regression}
\keyword{classification}
