#ifndef GRBPY_INTEGRATE
#define GRBPY_INTEGRATE

double trap(double (*f)(double, void *), double xa, double xb, int N, void *args);
double simp(double (*f)(double, void *), double xa, double xb, int N, void *args);
double romb(double (*f)(double, void *), double xa, double xb, int N, double atol, 
                double rtol, void *args);
void simp_v2(void (*f)(double, double *, double *, double *, int, void *),
                        double *I, double *t1, double *t2, int Nt, double xa, 
                        double xb, int N, void *args);

#endif
