#ifndef GRBPY_SHOCK
#define GRBPY_SHOCK

void Rudot2D(double t, double *x, void *argv, double *xdot);
void shockEvolveRK4(double *t, double *R, double *u, int N, double R0, 
                    double u0, void *args);

#endif
