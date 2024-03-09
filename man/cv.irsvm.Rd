\name{cv.irsvm}
\alias{cv.irsvm}
\alias{cv.irsvm.formula}
\alias{cv.irsvm.matrix}
\alias{cv.irsvm.default}
\title{Cross-validation for irsvm}
\description{Does k-fold cross-validation for irsvm}
\usage{
\method{cv.irsvm}{formula}(formula, data, weights, contrasts=NULL, ...)
\method{cv.irsvm}{matrix}(x, y, weights, ...)
\method{cv.irsvm}{default}(x,  ...)
}
\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{x}{\code{x} matrix as in \code{irsvm}.}
  \item{y}{response \code{y} as in \code{irsvm}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
  \item{contrasts}{the contrasts corresponding to \code{levels} from the respective models.}
  \item{\dots}{Other arguments that can be passed to \code{irsvm}.}
}
\details{Does a K-fold cross-validation to determine optimal tuning parameters in SVM: \code{cost} and \code{gamma} if \code{kernel} is nonlinear. It can also choose \code{s} used in \code{cfun}. 
}
\value{An object contains a list of ingredients of cross-validation including optimal tuning parameters.
  \item{residmat}{matrix with row values for \code{kernel="linear"} are \code{s, cost, error, k}, where \code{k} is the number of cross-validation fold. For nonlinear kernels, row values are \code{s, gamma, cost, error, k}.}
  \item{cost}{a value of \code{cost} that gives minimum cross-validated value in \code{irsvm}.} 
  \item{gamma}{a value of \code{gamma} that gives minimum cross-validated value in \code{irsvm}}
  \item{s}{value of \code{s} for \code{cfun} that gives minimum cross-validated value in \code{irsvm}.}
}

\examples{
\dontrun{
x <- matrix(rnorm(40*2), ncol=2)
y <- c(rep(-1, 20), rep(1, 20))
x[y==1,] <- x[y==1, ] + 1
irsvm.opt <- cv.irsvm(x, y, type="C-classification", s=1, kernel="linear", cfun="acave")
irsvm.opt$cost
irsvm.opt$gamma
irsvm.opt$s
}
}
\references{
 Zhu Wang (2024)
     \emph{Unified Robust Estimation}, \emph{Australian & New Zealand Journal         of      Statistics}. 66(1):77-102.
}
\author{Zhu Wang <zwang145@uthsc.edu>}

\seealso{\code{\link{irsvm}}}
\keyword{models}
\keyword{regression}

