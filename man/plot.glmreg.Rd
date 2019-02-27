\name{plot.glmreg}
\alias{plot.glmreg}
\title{plot coefficients from a "glmreg" object}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted \code{"glmreg"} object.
}
\usage{
\method{plot}{glmreg}(x, xvar = c("norm", "lambda", "dev"), label = FALSE, shade=TRUE, ...)
}
\arguments{
  \item{x}{fitted \code{"glmreg"} model}
  \item{xvar}{What is on the X-axis. \code{"norm"} plots against the
    L1-norm of the coefficients, \code{"lambda"} against the log-lambda
    sequence, and \code{"dev"} against the percent deviance explained.}
  \item{label}{If \code{TRUE}, label the curves with variable sequence
    numbers.}
  \item{shade}{Should nonconvex region be shaded?  Default is TRUE. Code developed for all \code{weights=1} only}
  \item{\dots}{Other graphical parameters to plot}
}
\details{
A coefficient profile plot is produced. 
}
\author{Zhu Wang wangz1@uthscsa.edu}
\seealso{\code{glmreg}, and \code{print}, \code{predict} and \code{coef} methods.}
\examples{
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
fit1=glmreg(x,y)
plot(fit1)
plot(fit1,xvar="lambda",label=TRUE)
}
\keyword{models}
\keyword{regression}
