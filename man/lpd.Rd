% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa-generics.R
\name{lpd}
\alias{lpd}
\title{Create generic \code{lpd} function for class}
\usage{
lpd(x)
}
\arguments{
\item{x}{An \code{Model} object.}
}
\value{
Numeric log-predictive density used in WAIC.
}
\description{
Create generic \code{lpd} function for class
}
\examples{
m1 <- Model(ndocs = 1)
lpd(m1)
}
