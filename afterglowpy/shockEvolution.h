#ifndef AFTERGLOWPY_SHOCK
#define AFTERGLOWPY_SHOCK

static const double T0_inj = 1.0e3;

double shockVel(double u);
double E_inj(double te, double L0, double q, double ts);
double L_inj(double te, double L0, double q, double ts);
void shockInitDecel(double t0, double *R0, double *u0, void *argv); 
void shockInitFind(double t0, double *R0, double *u0, double tRes, void *argv);
void Rudot2D(double t, double *x, void *argv, double *xdot);
void RuThdot3D(double t, double *x, void *argv, double *xdot, int spread);
void shockEvolveRK4(double *t, double *R, double *u, int N, double R0, 
                    double u0, void *args);
void shockEvolveSpreadRK4(double *t, double *R, double *u, double *th, int N, 
                            double R0, double u0, double th0, void *args,
                            int spread);

#endif
