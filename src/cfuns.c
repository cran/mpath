#include <R.h>
#include <Rmath.h>

int F77_SUB(cisnan)(double *x) { return ISNAN(*x); }
double F77_SUB(rlgamma)(double *x) { return lgammafn(*x); }
double F77_SUB(rdigamma)(double *x) { return digamma(*x); }
double F77_SUB(rtrigamma)(double *x) { return trigamma(*x); }
