\name{stan}
\alias{stan}

\title{standardize variables}
\description{
 Standardize variables. For each column, return mean 0 and mean value of sum of squares = 1.
  }
\usage{stan(x, weights)}

\arguments{
  \item{x}{numeric variables, can be a matrix or vector}
  \item{weights}{numeric positive vector of weights}
}
\value{A list with the following items.
\item{x}{standardized variables with each column: mean value 0 and mean value of sum of squares = 1.}
\item{meanx}{a vector of means for each column in the original \code{x}}
\item{normx}{a vector of scales for each column in the original \code{x}}
}
\author{Zhu Wang <wangz1@uthscsa.edu>}

