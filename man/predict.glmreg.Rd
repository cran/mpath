\name{predict.glmreg}
\alias{predict.glmreg}
\alias{coef.glmreg}
\title{Model predictions based on a fitted "glmreg" object.}
\description{
This function returns predictions from
a fitted \code{"glmreg"} object.
}
\usage{
\method{predict}{glmreg}(object,newx,newoffset,which=1:length(object$lambda),
type=c("link","response","class","coefficients","nonzero"), na.action=na.pass, ...)
\method{coef}{glmreg}(object,which=1:length(object$lambda),...)
}
\arguments{
  \item{object}{Fitted \code{"glmreg"} model object.}
  \item{newx}{Matrix of values at which predictions are to be made.  Not
    used for \code{type="coefficients"}}
  \item{which}{Indices of the penalty parameter \code{lambda} at which
    predictions are required.  By default, all indices are returned.}
  \item{type}{Type of prediction: \code{"link"} returns the
    linear predictors; \code{"response"} gives the fitted values;
    \code{"class"} returns the binomial outcome with the highest
    probability; \code{"coefficients"} returns the coefficients.}
  \item{newoffset}{an offset term used in prediction}
  \item{na.action}{action for missing data value}
  \item{\dots}{ arguments for predict }}
\value{The returned object depends on \code{type}.}
\references{
Zhu Wang, Shuangge Ma, Michael Zappitelli, Chirag Parikh, Ching-Yun Wang and Prasad Devarajan (2014)
\emph{Penalized Count Data Regression with Application to Hospital Stay after Pediatric Cardiac Surgery}, \emph{Statistical Methods in Medical Research}. 2014 Apr 17. [Epub ahead of print]
}
\author{Zhu Wang <wangz1@uthscsa.edu>}
\seealso{\code{\link{glmreg}}}
\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
fit <- glmreg(counts ~ outcome + treatment, data=d.AD, family="poisson")
predict(fit, newx=d.AD[,1:2])
summary(fit)
coef(fit)
}
\keyword{models}
\keyword{regression}
