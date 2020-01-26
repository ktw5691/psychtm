#ifndef COUNT_TOPIC_WORD_H
#define COUNT_TOPIC_WORD_H

#include <RcppArmadillo.h>

arma::mat count_topic_word(uint16_t K, uint32_t V, const arma::umat& doc_topic,
                           const arma::umat& doc_word);

#endif
