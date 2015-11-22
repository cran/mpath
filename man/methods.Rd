\name{methods}

\alias{AIC.glmreg}
\alias{AIC.zipath}
\alias{BIC.glmreg}
\alias{BIC.zipath}

\title{ Methods for mpath Objects }
\description{
  Methods for models fitted by coordinate descent algorithms.
}
\usage{
\method{AIC}{glmreg}(object, ..., k)
\method{AIC}{zipath}(object, ..., k)
\method{BIC}{glmreg}(object, ...)
\method{BIC}{zipath}(object, ...)

}
\arguments{
  \item{object}{ objects of class \code{glmreg} or \code{zipath}.}
 \item{\dots}{ additional arguments passed to callies. }
  \item{k}{  numeric, the \emph{penalty} per parameter to be used; the default
             \code{k = 2} is the classical AIC. \code{k} has been hard coded in the function and there is no impact to the value of \code{AIC} if \code{k} is changed }
}
\references{
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)      \emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]

Zhu Wang, Shuangge Ma, Ching-Yun Wang, Michael Zappitelli, Prasad Devarajan and Chirag R. Parikh (2014)
\emph{EM for Regularized Zero Inflated Regression Models with Applications to Postoperative Morbidity after Cardiac Surgery in Children}, \emph{Statistics in Medicine}. 33(29):5192-208.

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}
\author{Zhu Wang <zwang@connecticutchildrens.org>}

\keyword{methods}
