\name{pval.zipath}
\alias{pval.zipath}

\title{compute p-values from penalized zero-inflated model with multi-split data}

\description{
  compute p-values from penalized zero-inflated Poisson, negative binomial and geometric model with multi-split data
  }
\usage{
pval.zipath(formula, data, weights, subset, na.action, offset, standardize=TRUE,
            family = c("poisson", "negbin", "geometric"),
            penalty = c("enet", "mnet", "snet"), gamma.count = 3, 
            gamma.zero = 3, prop=0.5, trace=TRUE, B=10, ...)
}

\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights. If \code{standardize=TRUE}, weights are renormalized to weights/sum(weights). If \code{standardize=FALSE}, weights are kept as original input}
  \item{subset}{subset of data}
  \item{na.action}{how to deal with missing data}
  \item{offset}{Not implemented yet}
  \item{standardize}{logical value, should variables be standardized?}
  \item{family}{family to fit \code{zipath}}
  \item{penalty}{penalty considered as one of \code{enet, mnet, snet}.}
 \item{gamma.count}{The tuning parameter of the \code{snet} or \code{mnet} penalty for the count part of model.}
  \item{gamma.zero}{The tuning parameter of the \code{snet} or \code{mnet} penalty for the zero part of model.}  
  \item{prop}{proportion of data split, default is 50/50 split}
  \item{trace}{logical value, if TRUE, print detailed calculation results}
  \item{B}{number of repeated multi-split replications} 
  \item{...}{Other arguments passing to \code{glmreg_fit}}
}
\details{
  compute p-values from penalized zero-inflated Poisson, negative binomial and geometric model with multi-split data
  }
\value{
  \item{count.pval}{raw p-values in the count component}
  \item{zero.pval}{raw p-values in the zero component}
  \item{count.pval.q}{Q value for the count component}
  \item{zero.pval.q}{Q value for the zero component}
}
\references{
Nicolai Meinshausen, Lukas Meier and Peter Buehlmann (2013)
\emph{p-Values for High-Dimensional Regression}, \emph{Journal of the American Statistical Association}, 104(488), 1671--1681.

Zhu Wang, Shuangge Ma, Ching-Yun Wang, Michael Zappitelli, Prasad Devarajan and Chirag R. Parikh (2014)
\emph{EM for Regularized Zero Inflated Regression Models with Applications to Postoperative Morbidity after Cardiac Surgery in Children}, \emph{Statistics in Medicine}. 33(29):5192-208.

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}

\author{Zhu Wang <wangz1@uthscsa.edu>}

\keyword{models}
\keyword{regression}

 
