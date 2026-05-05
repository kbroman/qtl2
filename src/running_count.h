// running count, for calc_hotspots
#ifndef RUNNING_COUNT_H
#define RUNNING_COUNT_H

#include <Rcpp.h>

Rcpp::IntegerVector running_count(const Rcpp::NumericVector pos,
                                  const Rcpp::NumericVector result_pos,
                                  double window);

#endif // RUNNING_COUNT_H
