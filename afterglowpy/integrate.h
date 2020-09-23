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
double trap(double (*f)(double, void *), double xa, double xb, int N,
            void *args, int (*errf)(void *));
double simp(double (*f)(double, void *), double xa, double xb, int N,
            void *args, int (*errf)(void *));
double romb(double (*f)(double, void *), double xa, double xb, int N,
            double atol, double rtol, void *args, int *Neval, double *eps,
            int verbose, int (*errf)(void *), double *pfa, double *pfb);

double trap_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, struct Mesh3 *mout, int verbose,
                  int (*errf)(void *), double *pfa, double *pfb);
double simp_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, struct Mesh5 *mout, int verbose,
                  int (*errf)(void *), double *pfa, double *pfb);
double trapNL_adapt(double (*f)(double, void *), double xa, double xb,int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, struct Mesh5 *mout, int verbose,
                  int (*errf)(void *), double *pfa, double *pfb);
double hybrid_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, int verbose, int (*errf)(void *),
                  double *pfa, double *pfb);
double cadre_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, int verbose, int (*errf)(void *),
                  double *pfa, double *pfb);
double gk49_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, int verbose, int (*errf)(void *));
double gk715_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, int verbose, int (*errf)(void *));
double gk1021_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, int verbose, int (*errf)(void *));

/*
 * Internal functions for trap_adapt and simp_adapt.
 */

double m_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                 int (*processInterval)(double (*f)(double, void*), void *,
                                         Interval *, int (*errf)(void *)),
                 int (*splitInterval)(double (*f)(double, void *), void *,
                                Interval *, Interval *, Interval *,
                                int (*errf)(void *)),
                 double atol, double rtol, void *args, int *Neval,
                 double *eps, struct Mesh *mout, int verbose,
                 int (*errf)(void *));
double m3_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                 int (*initInterval)(double (*f)(double, void*), void *,
                                         Interval3 *, int (*errf)(void *),
                                     double *pfa, double *pfb),
                 int (*processInterval)(double (*f)(double, void*), void *,
                                         Interval3 *, int (*errf)(void *)),
                 int (*splitInterval)(double (*f)(double, void *), void *,
                                Interval3 *, Interval3 *, Interval3 *,
                                int (*errf)(void *)),
                 double atol, double rtol, void *args, int *Neval,
                 double *eps, struct Mesh3 *mout, int verbose,
                 int (*errf)(void *), double *pfa, double *pfb);
double m5_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                 int (*initInterval)(double (*f)(double, void*), void *,
                                         Interval5 *, int (*errf)(void *),
                                     double *pfa, double *pfb),
                 int (*processInterval)(double (*f)(double, void*), void *,
                                         Interval5 *, int (*errf)(void *)),
                 int (*splitInterval)(double (*f)(double, void *), void *,
                                Interval5 *, Interval5 *, Interval5 *,
                                int (*errf)(void *)),
                 double atol, double rtol, void *args, int *Neval,
                 double *eps, struct Mesh5 *mout, int verbose,
                 int (*errf)(void *), double *pfa, double *pbf);
double m9_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                 int (*initInterval)(double (*f)(double, void*), void *,
                                         Interval9 *, int (*errf)(void *),
                                     double *pfa, double *pfb),
                 int (*processInterval)(double (*f)(double, void*), void *,
                                         Interval9 *, int (*errf)(void *)),
                 int (*splitInterval)(double (*f)(double, void *), void *,
                                Interval9 *, Interval9 *, Interval9 *,
                                int (*errf)(void *)),
                 double atol, double rtol, void *args, int *Neval,
                 double *eps, struct Mesh9 *mout, int verbose,
                 int (*errf)(void *), double *pfa, double *pbf);

int trapInitInterval(double (*f)(double, void *), void *args, Interval3 *i,
                        int (*errf)(void *), double *pfa, double *pfb);
int trapProcessInterval(double (*f)(double, void *), void *args, Interval3 *i,
                        int (*errf)(void *));
int trapSplitInterval(double (*f)(double, void *), void *args,
                        Interval3 *i0, Interval3 *i1, Interval3 *i2,
                        int (*errf)(void *));

int simpInitInterval(double (*f)(double, void *), void *args, Interval5 *i,
                        int (*errf)(void *), double *pfa, double *pfb);
int simpProcessInterval(double (*f)(double, void *), void *args, Interval5 *i,
                        int (*errf)(void *));
int simpSplitInterval(double (*f)(double, void *), void *args,
                      Interval5 *i0, Interval5 *i1, Interval5 *i2,
                      int (*errf)(void *));

int trapNLInitInterval(double (*f)(double, void *), void *args, Interval5 *i,
                        int (*errf)(void *), double *pfa, double *pfb);
int trapNLProcessInterval(double (*f)(double, void *), void *args,
                          Interval5 *i, int (*errf)(void *));
int trapNLSplitInterval(double (*f)(double, void *), void *args,
                        Interval5 *i0, Interval5 *i1, Interval5 *i2,
                        int (*errf)(void *));

int cadreInitInterval(double (*f)(double, void *), void *args, Interval9 *i,
                        int (*errf)(void *), double *pfa, double *pfb);
int cadreProcessInterval(double (*f)(double, void *), void *args,
                          Interval9 *i, int (*errf)(void *));
int cadreSplitInterval(double (*f)(double, void *), void *args,
                        Interval9 *i0, Interval9 *i1, Interval9 *i2,
                        int (*errf)(void *));

int gk49ProcessInterval(double (*f)(double, void *), void *args,
                          Interval *i, int (*errf)(void *));
int gk49SplitInterval(double (*f)(double, void *), void *args,
                        Interval *i0, Interval *i1, Interval *i2,
                        int (*errf)(void *));

int gk715ProcessInterval(double (*f)(double, void *), void *args,
                          Interval *i, int (*errf)(void *));
int gk715SplitInterval(double (*f)(double, void *), void *args,
                        Interval *i0, Interval *i1, Interval *i2,
                        int (*errf)(void *));

int gk1021ProcessInterval(double (*f)(double, void *), void *args,
                          Interval *i, int (*errf)(void *));
int gk1021SplitInterval(double (*f)(double, void *), void *args,
                        Interval *i0, Interval *i1, Interval *i2,
                        int (*errf)(void *));

int gk_compute(double (*f)(double, void *), void *args, int (*errf)(void *),
               double c, double z0, const double xg[], const double xk[],
               const double wg[], const double wgk[], int ng,
               double *I, double *err);
#endif
