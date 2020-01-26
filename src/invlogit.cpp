#include <RcppArmadillo.h>

//' @title Compute inverse logit
//'
//' @name invlogit
//' @param x A double
//'
//' @return Inverse-logit(x) on probability scale.
double invlogit(double x) {
  double temp = exp(x) / (1.0 + exp(x));
  return temp;
}
