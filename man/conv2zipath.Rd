\name{conv2zipath}
\alias{conv2zipath}

\title{convert zeroinfl object to class zipath}
\description{convert zeroinfl object to class zipath, which then can be used to predict new data}

\usage{
conv2zipath(object, family=c("poisson", "negbin", "geometric"))
}
\arguments{
  \item{object}{an object of class zeroinfl}
  \item{family}{one of families in zeroinfl class}
}
\value{an object of class zipath}
\author{Zhu Wang <wangz1@uthscsa.edu>}

