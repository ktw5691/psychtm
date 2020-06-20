#ifndef GET_LOGPOST_H
#define GET_LOGPOST_H

#include <RcppArmadillo.h>

double get_lpost_eta(double ll, const arma::colvec& eta,
                     const arma::colvec& mu0, const arma::mat& sigma0);

double get_lpost_mlr(double ll,
                     const arma::colvec& eta, double sigma2,
                     const arma::colvec& mu0, const arma::mat& sigma0,
                     double a0, double b0);

double get_lpost_lda(double ll, const arma::mat& theta, const arma::mat& beta,
                     double gamma_, double alpha_);

double get_lpost_slda_norm(double ll, const arma::colvec& eta, double sigma2,
                           const arma::mat& theta, const arma::mat& beta,
                           const arma::colvec& mu0, const arma::mat& sigma0,
                           double gamma_, double alpha_,
                           double a0, double b0);

double get_lpost_slda_logit(double ll, const arma::colvec& eta,
                            const arma::mat& theta, const arma::mat& beta,
                            const arma::colvec& mu0, const arma::mat& sigma0,
                            double gamma_, double alpha_);

double eta_logpost_logit(const arma::mat& w, const arma::colvec& y,
                         const arma::vec& eta,
                         const arma::vec& mu0, const arma::mat& sigma0);

#endif
