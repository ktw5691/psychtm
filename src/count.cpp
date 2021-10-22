#include <RcppArmadillo.h>

//' @title Update number of times topic was drawn in document excluding current word
//'
//' @name count_topicd
//' @param topic The current topic index.
//' @param ndk A vector of the current number of draws of each topic in document d.
//'
//' @return A vector of the current number of draws of each topic in document d
//'   excluding word n.
//'
//' @noRd
void count_topicd(uint16_t topic, uint16_t doc, arma::mat& ndk) {
  // Exclude word n from topic counts in doc d
  ndk(doc, topic - 1)--;
  if (ndk(doc, topic - 1) < 0) ndk(doc, topic - 1) = 0;
}

//' @title Update number of times topic was drawn in corpus excluding current word
//'
//' @name count_topic_corpus
//' @param topic The current topic index.
//' @param nk A vector of the current number of draws of each topic in the corpus.
//'
//' @return A vector of the current number of draws of each topic in the corpus
//'   excluding the current word.
//'
//' @noRd
void count_topic_corpus(uint16_t topic, arma::vec& nk) {
  // Exclude word n from topic counts in corpus
  nk(topic - 1)--;
  if (nk(topic - 1) < 0) nk(topic - 1) = 0;
}

//' @title Update number of times word and topic co-occur in corpus
//'
//' @param word The current word index.
//' @param topic The current topic index.
//' @param nkm A K x V matrix of the current number of co-occurences of each topic
//'   and vocabulary term in the corpus.
//'
//' @return A K x V matrix of the current number of co-occurences of each topic
//'   and vocabulary term in the corpus excluding the current word.
//'
//' @noRd
void count_word_topic(uint32_t word, uint16_t topic, arma::mat& nkm) {
  // Exclude word n from topic-word counts
  nkm(topic - 1, word - 1)--;
  // Fix possible negative counts
  if (nkm(topic - 1, word - 1) < 0) nkm(topic - 1, word - 1) = 0;
}

//' @title Update all topic and topic-word counts in document and corpus
//'
//' @name update_zcounts
//' @param d The current document index.
//' @param word The current word index.
//' @param topic The current topic index.
//' @param ndk A vector of the current number of draws of each topic in document d.
//' @param nk A vector of the current number of draws of each topic in the corpus.
//' @param nkm A K x V matrix of the current number of co-occurences of each topic
//'   and vocabulary term in the corpus.
//'
//' @noRd
void update_zcounts(uint32_t d, uint32_t word, uint16_t topic, uint32_t doc,
                    arma::mat& ndk, arma::vec& nk, arma::mat& nkm) {

  count_topicd(topic, doc, ndk);
  count_topic_corpus(topic, nk);
  count_word_topic(word, topic, nkm);
}
