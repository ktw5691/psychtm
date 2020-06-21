#ifndef RDIRICHLET_CPP_H
#define RDIRICHLET_CPP_H

#include <RcppArmadillo.h>

arma::mat rdirichlet_cpp(uint32_t n, const arma::rowvec& alpha_);

#endif
