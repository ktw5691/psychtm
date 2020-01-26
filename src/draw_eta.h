#ifndef DRAW_ETA_H
#define DRAW_ETA_H

#include <RcppArmadillo.h>

arma::colvec draw_eta_norm(const arma::mat& w, const arma::vec& y,
                           long double sigma2, const arma::vec& mu0,
                           const arma::mat& sigma0);

arma::colvec draw_eta_logit(const arma::mat& w, const arma::colvec& y,
                            const arma::colvec& eta_prev,
                            const arma::colvec& mu0, const arma::mat& sigma0,
                            const arma::vec& proposal_sd,
                            arma::vec& attempt, arma::vec& accept);

#endif
