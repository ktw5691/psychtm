#ifndef DRAW_ETA_NORM_H
#define DRAW_ETA_NORM_H

#include <RcppArmadillo.h>

arma::colvec draw_eta_norm(const arma::mat& w, const arma::vec& y,
                           long double sigma2, const arma::vec& mu0,
                           const arma::mat& sigma0);

#endif
