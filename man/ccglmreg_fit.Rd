\name{ccglmreg_fit}
\alias{ccglmreg_fit}
\title{Internal function for penalized CC-estimators}
\description{
    Fit a linear model in penalized CC-family.  The
  regularization path is computed for the lasso (or elastic net penalty), scad (or snet) and mcp (or mnet penalty), at a grid
  of values for the regularization parameter lambda.
  }
\usage{
ccglmreg_fit(x, y, weights, offset, cfun="ccave", dfun="gaussian", s=NULL, 
             delta=0.1, fk=NULL, iter=10, reltol=1e-5, 
             penalty=c("enet","mnet","snet"), nlambda=100, lambda=NULL, 
             type.path=c("active", "nonactive"), decreasing=TRUE, 
             lambda.min.ratio=ifelse(nobs<nvars,.05, .001), alpha=1, gamma=3,
             rescale=TRUE, standardize=TRUE, intercept=TRUE, 
             penalty.factor= NULL, maxit=1000, type.init=c("bst", "co", "heu"), 
             init.family=NULL, mstop.init=10, nu.init=0.1, 
             eps=.Machine$double.eps, epscycle=10, thresh=1e-6, parallel=FALSE,
             n.cores=2, theta, trace=FALSE, tracelevel=1)
}

\arguments{
  \item{x}{input matrix, of dimension nobs x nvars; each row is an
  observation vector.} 
  \item{y}{response variable. Quantitative for \code{dfun=1} and -1/1 otherwise for classifications.}
  \item{weights}{observation weights. Can be total counts if responses are proportion matrices. Default is 1 for each observation}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of cases. Currently only one offset term can be included in the formula.}
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
  \item{dfun}{ character, type of convex downward function. \cr
    Valid options are:
    \itemize{
        \item \code{"gaussian"}
        \item \code{"gaussianC"}
        \item \code{"binomial"}
        }
    }
\item{s}{ tuning parameter of \code{cfun}. \code{s > 0} and can be equal    to 0 for \code{cfun="tcave"}. If \code{s} is too close to 0 for               \code{cfun="acave", "bcave", "ccave"}, the calculated weights can become 0    for all observations, thus crash the program.}
  \item{delta}{a small positive number provided by user only if                \code{cfun="gcave"} and \code{0 < s <1}}
  \item{fk}{ predicted values at an iteration in the COCO algorithm }
  \item{nlambda}{The number of \code{lambda} values - default is 100. The sequence may be truncated before \code{nlambda} is reached if a close to saturated model is fitted. See also \code{satu}.}
  \item{lambda}{ by default, the algorithm provides a sequence of regularization values, or a user supplied \code{lambda} sequence }
  \item{type.path}{solution path for \code{parallel=FALSE}. If \code{type.path="active"}, then cycle through only the active set in the next increasing \code{lambda} sequence. If \code{type.path="nonactive"}, no active set for each element of the lambda sequence and cycle through all the predictor variables.} 
  \item{lambda.min.ratio}{Smallest value for \code{lambda}, as a fraction of
    \code{lambda.max}, the (data derived) entry value (i.e. the smallest
  value for which all coefficients are zero except the intercept). Note, there is no closed formula for \code{lambda.max}. The default of \code{lambda.min.ratio} depends on the
  sample size \code{nobs} relative to the number of variables
  \code{nvars}. If \code{nobs > nvars}, the default is \code{0.001},
  close to zero.  If \code{nobs < nvars}, the default is \code{0.05}.}
  \item{alpha}{The \eqn{L_2} penalty mixing parameter, with
    \eqn{0 \le alpha\le 1}. \code{alpha=1} is lasso (mcp, scad) penalty; and \code{alpha=0} the ridge penalty. However, if \code{alpha=0}, one must provide \code{lambda} values.}
  \item{gamma}{The tuning parameter of the \code{snet} or \code{mnet} penalty.}
  \item{rescale}{logical value, if TRUE, adaptive rescaling of the penalty       parameter for \code{penalty="mnet"} or \code{penalty="snet"} with               \code{dfun="binomial"}. See \code{glmreg_fit}}
  \item{standardize}{logical value for x variable standardization, prior to
    fitting the model sequence. The coefficients are always returned on
    the original scale. Default is \code{standardize=TRUE}.}
  \item{intercept}{logical value: if TRUE (default), intercept(s) are fitted; otherwise, intercept(s)  are set to zero}
  \item{penalty.factor}{This is a number that multiplies \code{lambda} to allow
    differential shrinkage of coefficients. Can be 0 for some variables, which implies
    no shrinkage, and that variable is always included in the
    model. Default is same shrinkage for all variables.}
 \item{type.init}{a method to determine the initial values. If \code{type.init="ncl"}, an intercept-only model as initial parameter and run \code{ccglmreg} regularization path forward from lambda_max to lambda_min. If \code{type.init="heu"}, heuristic initial parameters and run \code{ccglmreg} path backward or forward depending on \code{decreasing}, between lambda_min and lambda_max. If \code{type.init="bst"}, run a boosting model with \code{bst} in package bst, depending on \code{mstop.init, nu.init} and run \code{ccglmreg} backward or forward depending on \code{decreasing}.}
 \item{init.family}{ character value for initial family, one of "clossR",     "closs","gloss","qloss", which can be used to derive an initial estimator, if   the selection is different from the default value }
 \item{mstop.init}{ an integer giving the number of boosting iterations when \code{type.init="bst"} }
 \item{nu.init}{ a small number (between 0 and 1) defining the step size or shrinkage parameter when \code{type.init="bst"}. }
 \item{decreasing}{ only used if \code{lambda=NULL}, a logical value used to determine regularization path direction either from lambda_max to a potentially modified lambda_min or vice versa if \code{type.init="bst", "heu"}. Since this is a nonconvex optimization, it is possible to generate different estimates for the same \code{lambda} depending on \code{decreasing}. The choice of \code{decreasing} picks different starting values.}
 \item{iter}{number of iteration in the COCO algorithm}
 \item{maxit}{Within each COCO algorithm iteration, maximum number of coordinate descent iterations for each \code{lambda} value; default is 1000.}
 \item{reltol}{convergency criteria in the COCO algorithm} 
 \item{eps}{If a coefficient is less than \code{eps} in magnitude, then it is reported to be 0}
 \item{epscycle}{ If \code{nlambda} > 1 and the relative loss values from two consecutive \code{lambda} values change > \code{epscycle}, then re-estimate parameters in an effort to avoid trap of local optimization.}
 \item{thresh}{Convergence threshold for coordinate descent. Defaults value is \code{1e-6}.}
 \item{penalty}{Type of regularization}
 \item{theta}{ an overdispersion scaling parameter for \code{family="negbin"}}
 \item{parallel, n.cores}{If \code{TRUE}, to compute solution of \code{lambda} with parallel computing in number of \code{n.cores}. If \code{FALSE}, sequential computing. If \code{NULL}, still sequential computing with a different convergence criteria based on penalized loss values}
 \item{trace, tracelevel}{If \code{TRUE}, fitting progress is reported. If \code{tracelevel=2}, deeper level of fitting progress is reported.}
}

\details{
A case weighted penalized least squares or logistic regression is fit by the COCO algorithm, where the loss function is a composite function of \code{cfun}o\code{type}.
  The sequence of robust models implied by \code{lambda} is fit by COCO along with coordinate
  descent. Note that the objective function is \deqn{weights*loss + \lambda*penalty,} if \code{standardize=FALSE} and \deqn{
  \frac{weights}{\sum(weights)}*loss + \lambda*penalty,} if \code{standardize=TRUE}. 
  }
\value{
An object with S3 class \code{"ccglmreg"} for the various types of models.
  \item{call}{the call that produced the model fit}
  \item{b0}{Intercept sequence of length \code{length(lambda)}}
  \item{beta}{A \code{nvars x length(lambda)} matrix of coefficients.}
  \item{lambda}{The actual sequence of \code{lambda} values used}
  \item{weights_update}{A \code{nobs x length(lambda)} matrix of weights computed by the COCO algorithm. The entry of i-th row and j-th column is the weight for the i-th observation and j-th \code{lambda} value.}
  \item{decreasing}{ if \code{lambda} is an increasing sequence or not, used to determine regularization path direction either from lambda_max to a potentially modified lambda_min or vice versa if \code{type.init="bst", "heu"}. }
}

\references{
   Zhu Wang (2020)
   \emph{Unified Robust Estimation via the COCO, arXiv e-prints}, \url{https://arxiv.org/abs/2010.02848}
  }

\author{Zhu Wang <wangz1@uthscsa.edu>}

\seealso{\code{\link{ccglmreg}}}
\keyword{models}
\keyword{regression}

 
