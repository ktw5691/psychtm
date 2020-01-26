#ifndef RMVNORM_CPP_H
#define RMVNORM_CPP_H

#include <RcppArmadillo.h>

arma::mat rmvnorm_cpp(uint32_t n, const arma::colvec& mu,
                      const arma::mat& sigma);
#endif
