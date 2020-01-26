#ifndef GET_LOG_NUMER_SAMPLEZ_H
#define GET_LOG_NUMER_SAMPLEZ_H

#include <RcppArmadillo.h>

arma::vec get_log_numer_samplez(uint32_t V, const arma::vec& ndk_n,
                                const arma::vec& nkm_n, const arma::vec& nk_n,
                                float alpha_, float gamma_);

#endif
