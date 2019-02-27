\name{rzi}
\alias{rzi}

\title{random number generation of zero-inflated count response}

\description{
random number generation of zero-inflated count response
}
\usage{
rzi(n, x, z, a, b, theta=1, family=c("poisson", "negbin", "geometric"), infl=TRUE)
}

\arguments{
  \item{n}{sample size of random number generation}
  \item{x}{design matrix of count model}
  \item{z}{design matrix of zero model}
  \item{a}{coefficient vector for \code{x}, length must be the same as column size of \code{x}}
  \item{b}{coefficient vector for \code{z}, length must be the same as column size of \code{z}}
  \item{theta}{dispersion parameter for \code{family="negbin"}}
  \item{family}{distribution of count model}
  \item{infl}{logical value, if TRUE, zero-inflated count response}
}
\details{random number generation of zero-inflated count response} 
\value{numeric vector of zero-inflated count response}
\references{
Zhu Wang, Shuangge Ma, Ching-Yun Wang, Michael Zappitelli, Prasad Devarajan and Chirag R. Parikh (2014)
\emph{EM for Regularized Zero Inflated Regression Models with Applications to Postoperative Morbidity after Cardiac Surgery in Children}, \emph{Statistics in Medicine}. 33(29):5192-208.

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}

\author{Zhu Wang <wangz1@uthscsa.edu>}

\keyword{models}
\keyword{regression}

 
