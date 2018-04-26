#ifndef GRBPY_SHOCK
#define GRBPY_SHOCK

double shockVel(double u);
void shockInitDecel(double t0, double *R0, double *u0, void *argv); 
void shockInitFind(double t0, double *R0, double *u0, double tRes, void *argv);
void Rudot2D(double t, double *x, void *argv, double *xdot);
void RuThdot3D(double t, double *x, void *argv, double *xdot);
void shockEvolveRK4(double *t, double *R, double *u, int N, double R0, 
                    double u0, void *args);
void shockEvolveSpreadRK4(double *t, double *R, double *u, double *th, int N, 
                            double R0, double u0, double th0, void *args);

#endif
