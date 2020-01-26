#ifndef GET_LOGLIKE_H
#define GET_LOGLIKE_H

#include <RcppArmadillo.h>

double get_ll_logit_yd(bool yd, double muhatd);
double get_ll_logit(const arma::colvec& y, const arma::mat& w,
                    const arma::colvec& eta);

double get_ll_mlr(const arma::colvec& y, const arma::mat& w,
                  const arma::colvec& eta, double sigma2);

double get_ll_lda(const arma::umat& zdocs, const arma::umat& docs,
                  const arma::mat& theta, const arma::mat& beta,
                  const Rcpp::IntegerVector& docs_index, const arma::colvec& N);

double get_ll_slda_norm(const arma::colvec& y, const arma::mat& w,
                        const arma::colvec& eta, double sigma2,
                        const arma::umat& zdocs, const arma::umat& docs,
                        const arma::mat& theta, const arma::mat& beta,
                        const Rcpp::IntegerVector& docs_index,
                        const arma::colvec& N);

double get_ll_slda_logit(const arma::colvec& y, const arma::mat& w,
                         const arma::colvec& eta,
                         const arma::umat& zdocs, const arma::umat& docs,
                         const arma::mat& theta, const arma::mat& beta,
                         const Rcpp::IntegerVector& docs_index,
                         const arma::colvec& N);

#endif
