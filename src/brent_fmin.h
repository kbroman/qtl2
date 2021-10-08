// 1-d optimization by Brent's method
#ifndef BRENT_FMIN_H
#define BRENT_FMIN_H

/***********************************************************************
 * Taken from R ver 3.2.2
 * R-3.2.2/src/library/stats/src/optimize.c
 ***********************************************************************/

// see brent_fmin.cpp

/* add attribute with clang to avoid clang-UBSAN warning on CRAN */
#if defined(__clang__)
#define QTL2_NO_SANITIZE __attribute__((no_sanitize("undefined")))
#else
#define QTL2_NO_SANITIZE
#endif

double QTL2_NO_SANITIZE qtl2_Brent_fmin(double ax, double bx, double (*f)(double, void *),
                                        void *info, double tol);

#endif // BRENT_FMIN_H
