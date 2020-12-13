\name{cv.ccsvm}
\alias{cv.ccsvm}
\alias{cv.ccsvm.formula}
\alias{cv.ccsvm.matrix}
\alias{cv.ccsvm.default}
\title{Cross-validation for ccsvm}
\description{Does k-fold cross-validation for ccsvm}
\usage{
\method{cv.ccsvm}{formula}(formula, data, weights, contrasts=NULL, ...)
\method{cv.ccsvm}{matrix}(x, y, weights, ...)
\method{cv.ccsvm}{default}(x,  ...)
}
\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{x}{\code{x} matrix as in \code{ccsvm}.}
  \item{y}{response \code{y} as in \code{ccsvm}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
  \item{contrasts}{the contrasts corresponding to \code{levels} from the respective models.}
  \item{\dots}{Other arguments that can be passed to \code{ccsvm}.}
}
\details{Does a K-fold cross-validation to determine optimal tuning parameters in SVM: \code{cost} and \code{gamma} if \code{kernel} is nonlinear. It can also choose \code{s} used in \code{cfun}. 
}
\value{An object contains a list of ingredients of cross-validation including optimal tuning parameters.
  \item{residmat}{matrix with row values for \code{kernel="linear"} are \code{s, cost, error, k}, where \code{k} is the number of cross-validation fold. For nonlinear kernels, row values are \code{s, gamma, cost, error, k}.}
  \item{cost}{a value of \code{cost} that gives minimum cross-validated value in \code{ccsvm}.} 
  \item{gamma}{a value of \code{gamma} that gives minimum cross-validated value in \code{ccsvm}}
  \item{s}{value of \code{s} for \code{cfun} that gives minimum cross-validated value in \code{ccsvm}.}
}

\examples{
\dontrun{
x <- matrix(rnorm(40*2), ncol=2)
y <- c(rep(-1, 20), rep(1, 20))
x[y==1,] <- x[y==1, ] + 1
ccsvm.opt <- cv.ccsvm(x, y, type="C-classification", s=1, kernel="linear", cfun="acave")
ccsvm.opt$cost
ccsvm.opt$gamma
ccsvm.opt$s
}
}
\references{
   Zhu Wang (2020) \emph{Unified Robust Estimation via the COCO, arXiv e-prints}, \url{https://arxiv.org/abs/2010.02848}
 }
\author{Zhu Wang <wangz1@uthscsa.edu>}

\seealso{\code{\link{ccsvm}}}
\keyword{models}
\keyword{regression}

