#ifndef EST_THETAD_H
#define EST_THETAD_H

#include <RcppArmadillo.h>

arma::rowvec est_thetad(const arma::rowvec& z_count, float alpha_);

#endif
