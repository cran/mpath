\name{ccsvm_fit}
\alias{ccsvm_fit}
\title{Fit iteratively re-weighted support vector machines for robust loss functions}
\description{
\code{ccsvm} is used to train a subject weighted support vector machine where the weights are provided iteratively from robust loss function with COCO algorithm. It can be used to carry out robust regression and binary classification.
}
\usage{
ccsvm_fit(x, y, weights, cfun="ccave", s=NULL, delta=0.0001, type = NULL, 
          kernel="radial", cost=1, epsilon = 0.1, iter=10, reltol=1e-5, 
          trace=FALSE, ...)
}
\arguments{
  \item{x}{a data matrix, a vector, or a sparse '\emph{design} matrix' (object of class
    \code{\link[Matrix]{Matrix}} provided by the \pkg{Matrix} package,
    or of class \code{\link[SparseM]{matrix.csr}}
    provided by the \pkg{SparseM} package, or of class
    \code{\link[slam]{simple_triplet_matrix}} provided by the \pkg{slam}
    package).}
  \item{y}{a response vector with one label for each row/component of
    \code{x}. Can be either a factor (for classification tasks)
    or a numeric vector (for regression).}
  \item{weights}{the weight of each subject. It should be in the same length of \code{y}.}
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
 \item{s}{ tuning parameter of \code{cfun}. \code{s > 0} and can be equal    to  0 for \code{cfun="tcave"}. If \code{s} is too close to 0 for \code{cfun="acave", "bcave", "ccave"}, the calculated weights can become 0 for all observations, thus crash the program.}
   \item{delta}{a small positive number provided by user only if                 \code{cfun="gcave"} and \code{0 < s <1}}
   \item{type}{\code{ccsvm} can be used as a classification
    machine, or as a regression machine.
    Depending of whether \code{y} is
    a factor or not, the default setting for \code{type} is \code{C-classification} or \code{eps-regression}, respectively, but may be overwritten by setting an explicit value.\cr
    Valid options are:
    \itemize{
      \item \code{C-classification}
      \item \code{nu-classification}
      \item \code{eps-regression}
      \item \code{nu-regression}
    }
  }
  \item{kernel}{the kernel used in training and predicting. You
     might consider changing some of the following parameters, depending
     on the kernel type.\cr
     \describe{
       \item{linear:}{\eqn{u'v}{u'*v}}
       \item{polynomial:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*v +      coef0)^degree}}
       \item{radial basis:}{\eqn{e^(-\gamma |u-v|^2)}{exp(-gamma*|u-v|^2)}}
       \item{sigmoid:}{\eqn{tanh(\gamma u'v + coef0)}{tanh(gamma*u'*v + coef0)}}
       }
     }
  \item{cost}{cost of constraints violation (default: 1)---it is the
     \sQuote{C}-constant of the regularization term in the Lagrange formulation. This is proportional to the inverse of \code{lambda} in \code{ccglmreg}.} 
  \item{epsilon}{epsilon in the insensitive-loss function (default: 0.1)}
\item{iter}{number of iteration in the COCO algorithm}
  \item{reltol}{convergency criteria in the COCO algorithm}
  \item{trace}{If \code{TRUE}, fitting progress is reported}
\item{\dots}{additional parameters for function \code{wsvm} in package \pkg{WeightSVM}}
}

\value{
  An object of class \code{"wsvm"} (see package \pkg{WeightSVM}) containing the fitted model, including:
  \item{SV}{The resulting support vectors (possibly scaled).}
  \item{index}{The index of the resulting support vectors in the data
    matrix. Note that this index refers to the preprocessed data (after
    the possible effect of \code{na.omit} and \code{subset})}
  \item{coefs}{The corresponding coefficients times the training labels.}
  \item{rho}{The negative intercept.}
  \item{sigma}{In case of a probabilistic regression model, the scale
    parameter of the hypothesized (zero-mean) laplace distribution estimated by
    maximum likelihood.}
  \item{probA, probB}{numeric vectors of length 2, number of
    classes, containing the parameters of the logistic distributions fitted to
    the decision values of the binary classifiers (1 / (1 + exp(a x + b))).}
}
\details{
    A case weighted SVM is fit by the COCO algorithm, where the loss function is a composite function of \code{cfun}o\code{type}, plus a \eqn{L\_2} penalty.
 Additional arguments include \code{degree, gamma, coef0}, 
 \code{class.weights, cachesize, tolerance, shrinking, propbability, fitted}, the same as \code{"wsvm"} in package \pkg{WeightSVM}.
}
\references{
   Zhu Wang (2020)
   \emph{Unified Robust Estimation via the COCO, arXiv e-prints}, \url{https://arxiv.org/abs/2010.02848}
  }
\author{
  Zhu Wang \email{wangz1@uthscsa.edu}
}
\seealso{\code{\link{print}}, \code{\link{predict}}, \code{\link{coef}} and     \code{\link{plot}}.}
\examples{
data(iris)
 iris <- subset(iris, Species \%in\% c("setosa", "versicolor"))
 # default with factor response:
  model <- ccsvm(Species ~ ., data = iris, kernel="linear", trace=TRUE)
  model <- ccsvm(Species ~ ., data = iris)
 # alternatively the traditional interface:
  x <- subset(iris, select = -Species)
  y <- iris$Species
model <- ccsvm(x, y)
  # test with train data
  pred <- predict(model, x)
  # (same as:)
  pred <- fitted(model)
 
  # Check accuracy:
  table(pred, y)
 # compute decision values and probabilities:
  pred <- predict(model, x, decision.values = TRUE)
  attr(pred, "decision.values")
 
  # visualize (classes by color, SV by crosses):
  plot(cmdscale(dist(iris[,-5])),
       col = as.integer(iris[,5]),
       pch = c("o","+")[1:100 \%in\% model$index + 1])
 
  ## try regression mode on two dimensions
 
  # create data
  x <- seq(0.1, 5, by = 0.05)
  y <- log(x) + rnorm(x, sd = 0.2)
 
  # estimate model and predict input values
  m   <- ccsvm(x, y)
  new <- predict(m, x)
 
  # visualize
 plot(x, y)
  points(x, log(x), col = 2)
  points(x, new, col = 4) 
}
 \keyword{neural}
 \keyword{nonlinear}
 \keyword{classif}
