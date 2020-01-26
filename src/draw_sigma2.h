#ifndef DRAW_SIGMA2_H
#define DRAW_SIGMA2_H

#include <RcppArmadillo.h>

long double draw_sigma2(float a0, float b0,
                        const arma::mat& w, const arma::colvec& y,
                        const arma::colvec& eta);

#endif
