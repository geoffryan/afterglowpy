#ifndef INTEGRATE_H
#define INTEGRATE_H

double trap(double (*f)(double, void *), double xa, double xb, int N, void *args);
double simp(double (*f)(double, void *), double xa, double xb, int N, void *args);
double romb(double (*f)(double, void *), double xa, double xb, int N, double atol, 
                double rtol, void *args);

#endif
