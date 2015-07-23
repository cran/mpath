\name{vuong.test}
\alias{vuong.test}
\title{Vuong's non-nested hypothesis test}
\description{
  Compares two models fit to the same data that do not nest via Vuong's
  non-nested test.
}

\usage{
vuong.test(m1, s1, m2, s2, type=c("None","AIC", "BIC"), digits = getOption("digits"))
}

\arguments{
  \item{m1}{model 1, an object inheriting from class \code{glm},
    \code{negbin}, \code{zeroinfl}, \code{glmreg}, \code{glmregNB} or \code{zipath}}
  \item{s1}{an integer value between 1 and m1$nlambda, only useful when model 1 is from 
\code{glmreg}, \code{glmregNB} or \code{zipath}}
  \item{m2}{model 2, as for model 1}
  \item{s2}{similar to s1}
  \item{type}{\code{type} of penalty on number of parameters when comparing models}
  \item{digits}{significant digits in printed result}
}
\details{
  The Vuong non-nested test is based on a comparison of the predicted
  probabilities of two models that do not nest, subject to penalty on number of parameters. The test can be used for penalized or non-penalized Poisson (negative binomial) regression models, zero-inflated Poisson (negative binomial) models.
}
\value{
  The test-statistic and \eqn{p} value
}
\references{Vuong, Q.H. 1989. Likelihood ratio tests for model
  selection and non-nested hypotheses. \emph{Econometrica}. \bold{57}, 307--333.

Achim Zeileis, Christian Kleiber and Simon Jackman (2008),                                                \emph{Regression Models for Count Data in R},
\emph{Journal of Statistical Software}, \bold{27(8)},
\url{http://www.jstatsoft.org/v27/i08/}
}
\author{Zhu Wang zwang@connecticutchildrens.org}

\examples{
\dontrun{
data("bioChemists", package="pscl")
## compare penalized Poisson GLM and ZIP
glm1 <- glmreg(art ~ ., data = bioChemists, x.keep=TRUE, family = "poisson")
zip <- zipath(art ~ . | ., data = bioChemists)
vuong.test(glm1, which.min(glm1$bic), zip, which.min(zip$bic))

## compare penalized negbin with zero-inflated negbin
nb1 <- glmregNB(art ~ ., data=bioChemists, x.keep=TRUE)
zinb <- zipath(art ~ . | ., data = bioChemists, family = "negbin")
vuong.test(nb1, which.min(nb1$bic), zinb, which.min(zinb$bic))
}
}

\keyword{models}
