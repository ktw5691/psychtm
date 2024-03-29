% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa-generics.R
\name{gamma_}
\alias{gamma_}
\title{Create generic \code{gamma_} function for class}
\usage{
gamma_(x)
}
\arguments{
\item{x}{An \code{Sldax} object.}
}
\value{
Double value of parameter for symmetric Dirichlet distribution prior
on the topic-word probabilities.
}
\description{
Create generic \code{gamma_} function for class
}
\examples{
m1 <- Sldax(ndocs = 1, nvocab = 2,
            topics = array(c(1, 2, 2, 1), dim = c(1, 4, 1)),
            theta = array(c(0.5, 0.5), dim = c(1, 2, 1)),
            beta = array(c(0.5, 0.5, 0.5, 0.5), dim = c(2, 2, 1)))
gamma_(m1)
}
