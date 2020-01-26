#ifndef TABLE_CPP_H
#define TABLE_CPP_H

#include <RcppArmadillo.h>

std::map<uint32_t, uint32_t> table_cpp(const arma::urowvec& v);

#endif
