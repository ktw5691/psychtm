% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa-generics.R
\name{eta<-}
\alias{eta<-}
\title{Create generic \verb{eta<-} function for class}
\usage{
eta(x) <- value
}
\arguments{
\item{x}{An \code{Model} object.}

\item{value}{Numeric vector of regression coefficients to assign to slot.}
}
\value{
None.
}
\description{
Create generic \verb{eta<-} function for class
}
\examples{
m1 <- Model(ndocs = 1)
eta(m1) <- matrix(c(-1.0, 1.0), nrow = 1, ncol = 2)
}
