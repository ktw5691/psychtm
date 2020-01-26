#ifndef GIBBS_SLDAX_CPP_H
#define GIBBS_SLDAX_CPP_H

#include <RcppArmadillo.h>

Rcpp::S4 gibbs_sldax_cpp(const arma::umat& docs, uint32_t V,
                         uint32_t m, uint32_t burn, uint32_t thin,
                         uint16_t K, uint8_t model,
                         const arma::colvec& y,
                         const arma::mat& x,
                         const arma::colvec& mu0,
                         const arma::mat& sigma0,
                         float a0, float b0,
                         arma::colvec eta_start,
                         arma::vec proposal_sd,
                         int interaction_xcol = -1,
                         float alpha_ = 0.1, float gamma_ = 1.01,
                         bool constrain_eta = true,
                         bool verbose = false, bool display_progress = false);

#endif
