#ifndef COUNT_H
#define COUNT_H

#include <RcppArmadillo.h>

void count_topicd(uint16_t topic, uint16_t doc, arma::mat& ndk);

void count_topic_corpus(uint16_t topic, arma::vec& nk);

void count_word_topic(uint32_t word, uint16_t topic, arma::mat& nkm);

void update_zcounts(uint32_t d, uint32_t word, uint16_t topic, uint32_t doc,
                    arma::mat& ndk, arma::vec& nk, arma::mat& nkm);

#endif
