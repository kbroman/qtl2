// running count, for calc_hotspots

#include "running_count.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".running_count")]]
IntegerVector running_count(const NumericVector pos,
                            const NumericVector result_pos,
                            double window)
{
    const int n_pos = pos.size();
    const int n_result = result_pos.size();
    IntegerVector result(n_result);

    window /= 2.0;

    for(int i=0; i<n_result; i++) {
        result[i] = 0;

        for(int j=0; j<n_pos; j++) {
            if(pos[j] >= result_pos[i]-window &&
               pos[j] <= result_pos[i]+window) result[i]++;
        }
    }

    return result;
}
