\name{gfunc}
\alias{gfunc}
\title{Convert response value to raw prediction in GLM}
\description{
 Compute response value to raw prediction such as linear predictor in GLM
}
\usage{
gfunc(mu, family, epsbino)
}
\arguments{
 \item{mu}{vector of numbers as response value in GLM, for instance, probability estimation if \code{family=2}}
 \item{family}{integer from 1-4, corresponding to "gaussian", "binomial", "poisson", "negbin", respectively}
 \item{epsbino}{a small positive value for \code{family=2} to avoid numeric unstability}
}

\value{linear predictor f=x'b for predictor x and coefficient b if the model is linear}

\keyword{regression}
\keyword{classification}
