\name{cv.ccsvm_fit}
\alias{cv.ccsvm_fit}
\title{Internal function of cross-validation for ccsvm}
\description{Internal function to conduct k-fold cross-validation for ccsvm}
\usage{
cv.ccsvm_fit(x, y, weights, cfun="ccave", s=c(1, 5), type=NULL, 
             kernel="radial", gamma=2^(-4:10), cost=2^(-4:4), 
             epsilon=0.1, balance=TRUE, nfolds=10, foldid, 
             trim_ratio=0.9, n.cores=2, ...)
}
\arguments{
 \item{x}{a data matrix, a vector, or a sparse '\emph{design} matrix' (object  of class
     \code{\link[Matrix]{Matrix}} provided by the \pkg{Matrix} package,
     or of class \code{\link[SparseM]{matrix.csr}}
     provided by the \pkg{SparseM} package, or of class
     \code{\link[slam]{simple_triplet_matrix}} provided by the \pkg{slam}
     package).}
   \item{y}{a response vector with one label for each row/component of
     \code{x}. Can be either a factor (for classification tasks)
     or a numeric vector (for regression).}
   \item{weights}{the weight of each subject. It should be in the same length    of \code{y}.}
    \item{cfun}{ character, type of convex cap (concave) function.\cr
    Valid options are:
     \itemize{
          \item \code{"hcave"}
         \item \code{"acave"}
         \item \code{"bcave"}
         \item \code{"ccave"}
         \item \code{"dcave"}
         \item \code{"ecave"}
         \item \code{"gcave"}
         \item \code{"tcave"}
      }
    }
\item{s}{ tuning parameter of \code{cfun}. \code{s > 0} and can be equal       to  0 for \code{cfun="tcave"}. If \code{s} is too close to 0 for                \code{cfun="acave", "bcave", "ccave"}, the calculated weights can become 0 for  all observations, thus crash the program.}
    \item{type}{\code{ccsvm} can be used as a classification
     machine, or as a regression machine.
     Depending of whether \code{y} is
     a factor or not, the default setting for \code{type} is \code{C-classification} or \code{eps-regression}, respectively, but may be overwritten  by setting an explicit value.\cr
     Valid options are:
     \itemize{
       \item \code{C-classification}
       \item \code{nu-classification}
       \item \code{eps-regression}
       \item \code{nu-regression}
     }
   }
   \item{kernel, gamma}{the kernel used in training and predicting. You
      might consider changing some of the following parameters, depending
      on the kernel type.\cr
      \describe{
        \item{linear:}{\eqn{u'v}{u'*v}}
        \item{polynomial:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*        v +      coef0)^degree}}
        \item{radial basis:}{\eqn{e^(-\gamma |u-v|^2)}{exp(-gamma*|u-v|^2)}}
        \item{sigmoid:}{\eqn{tanh(\gamma u'v + coef0)}{tanh(gamma*u'*v +         coef0)}}
      }
   }
   \item{cost}{cost of constraints violation (default: 1)---it is the
      \sQuote{C}-constant of the regularization term in the Lagrange formulation. This is proportional to the inverse of \code{lambda} in \code{ccglmreg}.}
   \item{epsilon}{epsilon in the insensitive-loss function (default: 0.1)}
   \item{balance}{for \code{type="C-classification", "nu-classification"} only}
  \item{nfolds}{number of folds >=3, default is 10}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in. If supplied,
    \code{nfold} can be missing and will be ignored.}
 \item{trim_ratio}{a number between 0 and 1 for trimmed least squares, useful if \code{type="eps-regression"} or \code{"nu-regression"}.} 
 \item{n.cores}{The number of CPU cores to use. The cross-validation loop
    will attempt to send different CV folds off to different cores.}
  \item{\dots}{Other arguments that can be passed to \code{ccsvm}.}
}
\details{
This function is the driving force behind \code{cv.ccsvm}. Does a K-fold cross-validation to determine optimal tuning parameters  in SVM: \code{cost} and \code{gamma} if \code{kernel} is nonlinear. It can also choose \code{s} used in \code{cfun}.    
}

\value{an object of class \code{"cv.ccsvm"} is returned, which is a
  list with the ingredients of the cross-validation fit.
\item{residmat}{matrix with row values for \code{kernel="linear"} are       \code{s, cost, error, k}, where \code{k} is the number of cross-validation fold. For nonlinear kernels, row values are \code{s, gamma, cost, error, k}. }
   \item{cost}{a value of \code{cost} that gives minimum cross-validated value in \code{ccsvm}.}
   \item{gamma}{a value of \code{gamma} that gives minimum cross-validated value in \code{ccsvm}}
   \item{s}{value of \code{s} for \code{cfun} that gives minimum cross-validated value in \code{ccsvm}.}
}
\references{
   Zhu Wang (2020)
   \emph{Unified Robust Estimation via the COCO, arXiv e-prints}, \url{https://arxiv.org/abs/2010.02848}
  }

\author{Zhu Wang <wangz1@uthscsa.edu>}
\seealso{\code{\link{cv.ccsvm}} and \code{\link{ccsvm}}}
\keyword{models}
\keyword{regression}

