\name{ncl_fit}
\alias{ncl_fit}
\title{Internal function to fit a nonconvex loss based robust linear model}
\description{
    Fit a linear model via penalized nonconvex loss function.  
  }
\usage{
ncl_fit(x,y, weights, cost=0.5, rfamily=c("clossR", "closs", "gloss", "qloss"), s=NULL, 
fk=NULL, iter=10, del=1e-10, trace=FALSE)
}

\arguments{
  \item{x}{input matrix, of dimension nobs x nvars; each row is an
  observation vector.} 
  \item{y}{response variable. Quantitative for \code{rfamily="clossR"} and -1/1 for classifications.}
  \item{weights}{observation weights. Can be total counts if responses are proportion matrices. Default is 1 for each observation}
  \item{cost}{ price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{rfamily}{Response type and relevant loss functions (see above)}
  \item{s}{ nonconvex loss tuning parameter for robust regression and classification. }
  \item{fk}{ predicted values at an iteration in the MM algorithm }
 \item{iter}{number of iteration in the MM algorithm}
 \item{del}{convergency critera} 
 \item{trace}{If \code{TRUE}, fitting progress is reported}
}
\details{
  The robust linear model is fit by majorization-minimization along with least squares.
  Note that the objective function is \deqn{1/2*weights*loss}.
  }
\value{
An object with S3 class \code{"ncl"} for the various types of models.
  \item{call}{the call that produced the model fit}
  \item{fitted.values}{predicted values}  
\item{h}{ pseudo response values in the MM algorithm }
}

\author{Zhu Wang <zwang@connecticutchildrens.org>}

\seealso{\code{\link{ncl}}}
\keyword{models}
\keyword{regression}

 
