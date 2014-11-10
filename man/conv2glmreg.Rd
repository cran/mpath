\name{conv2glmreg}
\alias{conv2glmreg}

\title{convert glm object to class glmreg}
\description{convert glm object to class glmreg, which then can be used for other purposes}

\usage{
conv2glmreg(object, family=c("poisson", "negbin"))
}
\arguments{
  \item{object}{an object of class glm}
  \item{family}{one of families in glm class}
}
\value{an object of class glmreg}
\author{Zhu Wang <zwang@connecticutchildrens.org>}

