\name{be.zeroinfl}
\alias{be.zeroinfl}

\title{conduct backward stepwise variable elimination for zero inflated count regression}

\description{
 conduct backward stepwise variable elimination for zero inflated count regression from zeroinfl function}

\usage{
be.zeroinfl(object, data, dist=c("poisson", "negbin", "geometric"), alpha=0.05, 
            trace=FALSE)
}

\arguments{
  \item{object}{an object from function zeroinfl}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{dist}{one of the distributions in zeroinfl function}
  \item{alpha}{significance level of variable elimination}
  \item{trace}{logical value, if TRUE, print detailed calculation results}
}
\details{
 conduct backward stepwise variable elimination for zero inflated count regression from zeroinfl function}
\value{an object of zeroinfl with all variables having p-values less than the significance level \code{alpha}}
\references{
Zhu Wang, Shuangge Ma, Ching-Yun Wang, Michael Zappitelli, Prasad Devarajan and Chirag R. Parikh (2014)
\emph{EM for Regularized Zero Inflated Regression Models with Applications to Postoperative Morbidity after Cardiac Surgery in Children}, \emph{Statistics in Medicine}. 33(29):5192-208.

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}

\author{Zhu Wang <wangz1@uthscsa.edu>}

\keyword{models}
\keyword{regression}

 
