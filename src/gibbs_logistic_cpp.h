#ifndef GIBBS_LOGISTIC_CPP_H
#define GIBBS_LOGISTIC_CPP_H

#include <RcppArmadillo.h>

Rcpp::S4 gibbs_logistic_cpp(uint32_t m, uint32_t burn, uint32_t thin,
                            const arma::colvec& y, const arma::mat& x,
                            const arma::colvec& mu0, const arma::mat& sigma0,
                            arma::colvec eta_start, arma::vec proposal_sd,
                            bool verbose = false, bool display_progress = false);

#endif
