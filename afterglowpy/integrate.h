#ifndef AFTERGLOWPY_INTEGRATE
#define AFTERGLOWPY_INTEGRATE

#include "interval.h"

/*
 * Various routines for integrating 1D functions.  
 * trap() and simp() are fixed stencil implementations of the Trapezoid Rule
 * and Simpson's Rule respectively.  
 * 
 * romb() is an adaptive Romberg integrator.
 * 
 * trap_adapt() and simp_adapt() are globally adaptive integrators based on 
 *    the Trapezoid and Simpson's rule, respectively. They successively bisect
 *    the integration domain into subintervals, prioritizing the subintervals
 *    with largest (estimated) error, until the total absolute error estimate
 *    is within tolerance.
 */

/*
 * Integration routines
 */
double trap(double (*f)(double, void *), double xa, double xb, int N, void *args);
double simp(double (*f)(double, void *), double xa, double xb, int N, void *args);
double romb(double (*f)(double, void *), double xa, double xb, int N,
            double atol, double rtol, void *args, int *Neval, double *eps,
            int verbose);

double trap_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, struct Mesh3 *mout, int verbose);
double simp_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, struct Mesh5 *mout, int verbose);

/*
 * Internal functions for trap_adapt and simp_adapt.
 */

double m3_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                 int (*processInterval)(double (*f)(double, void*), void *,
                                         Interval3 *),
                 int (*splitInterval)(double (*f)(double, void *), void *,
                                Interval3 *, Interval3 *, Interval3 *),
                 double atol, double rtol, void *args, int *Neval,
                 double *eps, struct Mesh3 *mout, int verbose);
double m5_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                 int (*processInterval)(double (*f)(double, void*), void *,
                                         Interval5 *),
                 int (*splitInterval)(double (*f)(double, void *), void *,
                                Interval5 *, Interval5 *, Interval5 *),
                 double atol, double rtol, void *args, int *Neval,
                 double *eps, struct Mesh5 *mout, int verbose);

int trapProcessInterval(double (*f)(double, void *), void *args, Interval3 *i);
int trapSplitInterval(double (*f)(double, void *), void *args,
                        Interval3 *i0, Interval3 *i1, Interval3 *i2);

int simpProcessInterval(double (*f)(double, void *), void *args, Interval5 *i);
int simpSplitInterval(double (*f)(double, void *), void *args,
                      Interval5 *i0, Interval5 *i1, Interval5 *i2);
#endif
