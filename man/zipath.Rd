\name{zipath}
\alias{zipath}
\alias{print.zipath}

\title{Fit zero-inflated count data linear model with lasso (or elastic net), snet or mnet regularization}
\description{
  Fit zero-inflated regression models for count data via penalized maximum likelihood.
}

\usage{
zipath(formula, data, weights, subset, na.action, offset,
standardize = TRUE, family = c("poisson", "negbin","geometric"), 
link = c("logit", "probit", "cloglog", "cauchit", "log"), 
penalty = c("enet", "mnet", "snet"), start = NULL, model = TRUE, 
y = TRUE, x = FALSE, nlambda = 100, lambda.count = NULL,lambda.zero = NULL, 
 penalty.factor.count=NULL, penalty.factor.zero=NULL,
lambda.count.min.ratio = .0001, lambda.zero.min.ratio = .1, 
alpha.count = 1, alpha.zero = alpha.count, gamma.count = 3, 
gamma.zero = gamma.count, rescale=FALSE, init.theta, theta.fixed=FALSE,
EM = TRUE, maxit.em=200, convtype=c("count", "both"), maxit = 1000, 
maxit.theta = 1, reltol = 1e-5, eps.bino=1e-5, shortlist=FALSE, trace = FALSE, ...)
}
\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{weights}{optional numeric vector of weights.}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor of the count model or zero model. See below for an example.}
\item{standardize}{Logical flag for x variable standardization, prior to
    fitting the model sequence. The coefficients are always returned on
    the original scale. Default is \code{standardize=TRUE}.}
  \item{family}{character specification of count model family (a log link is 
    always used).}
  \item{link}{character specification of link function in the binary
    zero-inflation model (a binomial family is always used).}
  \item{model, y, x}{logicals. If \code{TRUE} the corresponding components
    of the fit (model frame, response, model matrix) are returned.}
  \item{penalty}{penalty considered as one of \code{enet, mnet, snet}.}
  \item{start}{starting values for the parameters in the linear predictor.}
  \item{nlambda}{number of \code{lambda} value, default value is 100. The sequence may be truncated before \code{nlambda} is reached if a close to saturated model for the zero component is fitted.}
  \item{lambda.count}{A user supplied \code{lambda.count} sequence. Typical usage
    is to have the
    program compute its own \code{lambda.count} and \code{lambda.zero} sequence based on
    \code{nlambda} and \code{lambda.min.ratio}.} 
  \item{lambda.zero}{A user supplied \code{lambda.zero} sequence.}
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
  \item{init.theta}{The initial value of \code{theta} for \code{family="negbin"}.} 
  \item{theta.fixed}{Logical value only used for \code{family="negbin"}. If \code{TRUE}, \code{theta} is not updated.}
\item{EM}{Using \code{EM} algorithm. Not implemented otherwise}
\item{convtype}{convergency type, default is for count component only for speedy computation} 
\item{maxit.em}{Maximum number of EM algorithm}
\item{maxit}{Maximum number of coordinate descent algorithm}
 \item{maxit.theta}{Maximum number of iterations for estimating \code{theta} scaling parameter if family="negbin". Default value \code{maxit.theta} may be increased, yet may slow the algorithm}
 \item{eps.bino}{a lower bound of probabilities to be claimed as zero, for computing weights and related values when \code{family="binomial"}.}
\item{reltol}{Convergence criteria, default value 1e-5 may be reduced to make more accurate yet slow}
\item{shortlist}{logical value, if TRUE, limited results return}
\item{trace}{If \code{TRUE}, progress of algorithm is reported}
\item{...}{Other arguments which can be passed to from \code{glmreg}}
}

\details{
  The algorithm fits penalized zero-inflated count data regression models using the coordinate descent algorithm within the EM algorithm.
  The returned fitted model object is of class \code{"zipath"} and is similar
  to fitted \code{"glm"} and \code{"zeroinfl"} objects. For elements such as \code{"coefficients"} a list is returned with elements for the zero and count component,
  respectively. For details see below.
  
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"zipath"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{predict}}. See \code{\link{predict.zipath}} for more details
  on all methods.

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
  \item{contrasts}{a list with elements \code{"count"} and \code{"zero"}
    containing the contrasts corresponding to \code{levels} from the
    respective models,}
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

\seealso{\code{\link[stats]{glm}},
  \code{\link{glmreg}}, \code{\link{glmregNB}}
}

\examples{
\dontrun{
## data
data("bioChemists", package = "pscl")

## without inflation
## ("art ~ ." is "art ~ fem + mar + kid5 + phd + ment")
fm_pois <- glmreg(art ~ ., data = bioChemists, family = "poisson")
coef(fm_pois)
fm_nb <- glmregNB(art ~ ., data = bioChemists)
coef(fm_nb)
## with simple inflation (no regressors for zero component)
fm_zip <- zipath(art ~ . | 1, data = bioChemists, nlambda=10)
summary(fm_zip)
fm_zinb <- zipath(art ~ . | 1, data = bioChemists, family = "negbin", nlambda=10)
summary(fm_zinb)
## inflation with regressors
## ("art ~ . | ." is "art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment")
fm_zip2 <- zipath(art ~ . | ., data = bioChemists, nlambda=10)
summary(fm_zip2)
fm_zinb2 <- zipath(art ~ . | ., data = bioChemists, family = "negbin", nlambda=10)
summary(fm_zinb2)
### non-penalized regression, compare with zeroinfl
fm_zinb3 <- zipath(art ~ . | ., data = bioChemists, family = "negbin", 
lambda.count=0, lambda.zero=0, reltol=1e-12)
summary(fm_zinb3)
fm_zinb4 <- zerofinfl(art ~ . | ., data = bioChemists, dist = "negbin")
summary(fm_zinb4)
### offset
exposure <- rep(0.5, dim(bioChemists)[1])
fm_zinb <- zipath(art ~ . +offset(log(exposure))| ., data = bioChemists, 
		  family = "poisson", nlambda=10)
print(coef(fm_zinb))
print(predict(fm_zinb))

}
}

\keyword{regression}
