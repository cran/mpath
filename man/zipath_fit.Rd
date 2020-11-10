\name{zipath_fit}
\alias{zipath_fit}

\title{Internal function to fit zero-inflated count data linear model with lasso (or elastic net), snet or mnet regularization}
\description{
  Fit zero-inflated regression models for count data via penalized maximum likelihood.
}

\usage{
zipath_fit(X, Z, Y, weights, offsetx, offsetz, standardize=TRUE, 
           intercept = TRUE, family = c("poisson", "negbin", "geometric"), 
           link = c("logit", "probit", "cloglog", "cauchit", "log"), 
           penalty = c("enet", "mnet",  "snet"), start = NULL, y = TRUE, 
           x = FALSE, nlambda=100, lambda.count=NULL, lambda.zero=NULL, 
           type.path=c("active", "nonactive"), penalty.factor.count=NULL, 
           penalty.factor.zero=NULL, lambda.count.min.ratio=.0001, 
           lambda.zero.min.ratio=.1, alpha.count=1, alpha.zero=alpha.count, 
           gamma.count=3, gamma.zero=gamma.count, rescale=FALSE, 
           init.theta=NULL, theta.fixed=FALSE, EM=TRUE, maxit.em=200, 
           convtype=c("count", "both"), maxit= 1000, maxit.theta =10, 
           reltol = 1e-5, thresh=1e-6, eps.bino=1e-5, shortlist=FALSE, 
           trace=FALSE, ...)
}
\arguments{
  \item{X}{predictor matrix of the count model}
  \item{Z}{predictor matrix of the zero model}
  \item{Y}{response variable}
  \item{weights}{optional numeric vector of weights.}
  \item{offsetx}{optional numeric vector with an a priori known component to be
    included in the linear predictor of the count model.}
  \item{offsetz}{optional numeric vector with an a priori known component to be
    included in the linear predictor of the zero model.}
  \item{intercept}{Should intercept(s) be fitted (default=TRUE) or set to
   zero (FALSE)}
  \item{standardize}{Logical flag for x variable standardization, prior to
    fitting the model sequence. The coefficients are always returned on
    the original scale. 
    Default is \code{standardize=TRUE}.}
  \item{family}{character specification of count model family (a log link is 
    always used).}
  \item{link}{character specification of link function in the binary
    zero-inflation model (a binomial family is always used).}
  \item{y, x}{logicals. If \code{TRUE} the corresponding response and model matrix are returned.}
  \item{penalty}{penalty considered as one of \code{enet, mnet, snet}.}
  \item{start}{starting values for the parameters in the linear predictor.}
  \item{nlambda}{number of \code{lambda} value, default value is 100. The sequence may be truncated before \code{nlambda} is reached if a close to saturated model for the zero component is fitted.}
  \item{lambda.count}{A user supplied \code{lambda.count} sequence. Typical usage is to have the program compute its own \code{lambda.count} and \code{lambda.zero} sequence based on
    \code{nlambda} and \code{lambda.min.ratio}.} 
  \item{lambda.zero}{A user supplied \code{lambda.zero} sequence.}
\item{type.path}{solution path with default value \code{"active"}, which is less time computing than \code{"nonactive"}. If \code{type.path="nonactive"}, no active set for each element of the lambda sequence and cycle through all the predictor variables. 
    If \code{type.path="active"}, then cycle through only the active set, then cycle through all the variables for the same penalty parameter. See details below.} 
  \item{penalty.factor.count, penalty.factor.zero}{These are numeric vectors with the same length as predictor variables. that multiply \code{lambda.count, lambda.zero}, respectively, to allow differential shrinkage of coefficients. Can be 0 for some variables, which implies
    no shrinkage, and that variable is always included in the
    model. Default is same shrinkage for all variables.}
 \item{lambda.count.min.ratio, lambda.zero.min.ratio}{Smallest value for \code{lambda.count} 
       and \code{lambda.zero}, respectively, as a fraction of
    \code{lambda.max}, the (data derived) entry value (i.e. the smallest
  value for which all coefficients are zero except the intercepts).
Note, there is a closed formula for \code{lambda.max} for \code{penalty="enet"}. If \code{rescale=TRUE}, \code{lambda.max} is the same for \code{penalty="mnet" or "snet"}. Otherwise, some modifications are required. In the current implementation, for small \code{gamma} value, the square root of the computed \code{lambda.zero[1]} is used when \code{penalty="mnet" or "snet"}. 
}
  \item{alpha.count}{The elastic net mixing parameter for the count part of model.}  
  \item{alpha.zero}{The elastic net mixing parameter for the zero part of model.}  
  \item{gamma.count}{The tuning parameter of the \code{snet} or \code{mnet} penalty for the count part of model.}
  \item{gamma.zero}{The tuning parameter of the \code{snet} or \code{mnet} penalty for the zero part of model.}
  \item{rescale}{logical value, if TRUE, adaptive rescaling}
  \item{init.theta}{The initial value of \code{theta} for \code{family="negbin"}. This is set to \code{NULL} since version 0.3-24.} 
  \item{theta.fixed}{Logical value only used for \code{family="negbin"}. If \code{TRUE} and \code{init.theta} is provided with a numeric value > 0, then \code{init.theta} is not updated. 
      If \code{theta.fixed=FALSE}, then \code{init.theta} will be updated. In this case, if \code{init.theta=NULL}, its initial value is computed with intercept-only zero-inflated negbin model. }
\item{EM}{Using \code{EM} algorithm. Not implemented otherwise}
\item{convtype}{convergency type, default is for count component only for speedy computation} 
\item{maxit.em}{Maximum number of EM algorithm}
\item{maxit}{Maximum number of coordinate descent algorithm}
 \item{maxit.theta}{Maximum number of iterations for estimating \code{theta} scaling parameter if family="negbin". Default value \code{maxit.theta} may be increased, yet may slow the algorithm}
 \item{eps.bino}{a lower bound of probabilities to be claimed as zero, for computing weights and related values when \code{family="binomial"}.}
\item{reltol}{Convergence criteria, default value 1e-5 may be reduced to make more accurate yet slow}
\item{thresh}{Convergence threshold for coordinate descent. Defaults value is \code{1e-6}.}
\item{shortlist}{logical value, if TRUE, limited results return}
\item{trace}{If \code{TRUE}, progress of algorithm is reported}
\item{...}{Other arguments which can be passed to \code{glmreg} or \code{glmregNB}}
}

\details{
  The algorithm fits penalized zero-inflated count data regression models using the coordinate descent algorithm within the EM algorithm.
  The returned fitted model object is of class \code{"zipath"} and is similar
  to fitted \code{"glm"} and \code{"zeroinfl"} objects. For elements such as \code{"coefficients"} a list is returned with elements for the zero and count component,
  respectively. 
 
  If \code{type.path="active"}, the algorithm iterates for a pair (lambda_count, lambda_zero) in a loop:\cr
  Step 1: For initial coefficients start_count of the count model and start_zero of the zero model, the EM algorithm is iterated until convergence for the active set with non-zero coefficients determined from start_count and start_zero, respectively.\cr
  Step 2: EM is iterated for all the predict variables once.\cr
  Step 3: If active set obtained from Step 2 is the same as in Step 1, stop; otherwise, repeat Step 1 and Step 2.\cr
  If \code{type.path="nonactive"}, the EM algorithm iterates for a pair (lambda_count, lambda_zero) with all the predict variables until convergence. 

  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"zipath"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{predict}}. 
  See \code{\link{predict.zipath}} for more details on all methods.

  The program may terminate with the following message:

\code{Error in: while (j <= maxit.em && !converged)} \code{\{ :}\cr
\code{ Missing value, where TRUE/FALSE is necessary}\cr
\code{Calls: zipath}\cr
\code{Additionally: Warning:}\cr
\code{In glmreg_fit(Znew, probi, weights = weights, standardize = standardize,  :}\cr
\code{ saturated model, exiting ...}\cr
\code{Execution halted}\cr

  One possible reason is that the fitted model is too complex for the data. There are two suggestions to overcome the error. One is to reduce the number of variables. Second, find out what lambda values caused the problem and omit them. Try with other lambda values instead.

}

\value{
  An object of class \code{"zipath"}, i.e., a list with components including
  \item{coefficients}{a list with elements \code{"count"} and \code{"zero"}
    containing the coefficients from the respective models,}
  \item{residuals}{a vector of raw residuals (observed - fitted),}
  \item{fitted.values}{a vector of fitted means,}
  \item{weights}{the case weights used,}
  \item{terms}{a list with elements \code{"count"}, \code{"zero"} and
    \code{"full"} containing the terms objects for the respective models,}
  \item{theta}{estimate of the additional \eqn{\theta}{theta} parameter of the
    negative binomial model (if a negative binomial regression is used),}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{family}{character string describing the count distribution used,}
  \item{link}{character string describing the link of the zero-inflation model,}
  \item{linkinv}{the inverse link function corresponding to \code{link},}
  \item{converged}{logical value, TRUE indicating successful convergence of \code{zipath}, FALSE indicating otherwise}
  \item{call}{the original function call}
  \item{formula}{the original formula}
  \item{levels}{levels of the categorical regressors}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the response count vector (if \code{y = TRUE}),}
  \item{x}{a list with elements \code{"count"} and \code{"zero"}
    containing the model matrices from the respective models
    (if \code{x = TRUE}),}
}

\references{     
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)      \emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]

Zhu Wang, Shuangge Ma, Ching-Yun Wang, Michael Zappitelli, Prasad Devarajan and Chirag R. Parikh (2014)
\emph{EM for Regularized Zero Inflated Regression Models with Applications to Postoperative Morbidity after Cardiac Surgery in Children}, \emph{Statistics in Medicine}. 33(29):5192-208.

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}

\author{Zhu Wang <wangz1@uthscsa.edu>}

\seealso{\code{\link{zipath}},
  \code{\link{glmreg}}, \code{\link{glmregNB}}
}

\keyword{regression}
