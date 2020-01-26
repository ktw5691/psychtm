#ifndef DRAW_ZDN_H
#define DRAW_ZDN_H

uint16_t draw_zdn(arma::vec& log_num);

uint16_t draw_zdn_lda(uint32_t V,
                      const arma::vec& ndk_n, const arma::vec& nkm_n,
                      const arma::vec& nk_n, float alpha_, float gamma_);

uint16_t draw_zdn_slda_norm(double yd, const arma::rowvec& w_d,
                            const arma::vec& eta, double sigma2,
                            uint32_t V, const arma::vec& ndk_n,
                            const arma::vec& nkm_n, const arma::vec& nk_n,
                            float alpha_, float gamma_);

uint16_t draw_zdn_slda_logit(double yd, const arma::rowvec& w_d,
                             const arma::vec& eta,
                             uint32_t V, const arma::vec& ndk_n,
                             const arma::vec& nkm_n, const arma::vec& nk_n,
                             float alpha_, float gamma_);

#endif
