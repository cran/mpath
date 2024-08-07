\name{mpath-internal}
\title{Internal mpath functions}
\alias{balanced.folds}
\alias{cv.folds}
\alias{findlam}
\alias{g}
\alias{init}
\alias{llfun}
\alias{pen_sum}
\alias{predictzeroinfl1}
\alias{cfun2num}
\alias{y2num}
\alias{y2num4glm}
\alias{check_s}
\description{Internal mpath functions}
\author{Zhu Wang, author of ncvreg, author of pamr:balanced.folds (see COPYRIGHTS)}
\details{These are not intended for use by users. \code{findlam} is to compute the maximum value of \code{lambda} so that the coefficient estimation is zero. \code{g} is to compute a link function (not used at this moment), \code{init} is to compute eta, which is the estimated beta_0 in the intercept-only model. \code{pen_sum} is to compute the total penalty value from the estimated coefficients, \code{llfun} to compute the log-likelihood function. \code{convexMin} (from ncvreg) calculates the
    lowest index for which the penalized objective function is locally convex. 

Functions not used anymore: pen_eval, pen_sum, pen_zipath and zipath.control
}

\keyword{internal}
