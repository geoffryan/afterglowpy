#ifndef AFTERGLOWPY_SHOCK
#define AFTERGLOWPY_SHOCK

enum{ENV_ISM, ENV_WIND, ENV_PL, ENV_STEP};

static const double T0_inj = 1.0e3;

struct shockParams
{
    int spread;
    double E0;
    double Mej;
    double L0_inj;
    double q_inj;
    double t0_inj;
    double ts_inj;
    double E0_refresh;
    double k_refresh;
    double umin_refresh;
    double thetaCore;
    double theta0;
    double thetaCoreGlobal;
    int envType;
    double rho0_env;
    double R0_env;
    double k_env;
    double rho1_env;
};

void setup_shockParams(struct shockParams *pars, int spread,
                       double E0, double Mej, 
                       int envType, double env_rho0, double R0_env,
                       double k_env, double rho1_env,
                       double L0_inj, double q_inj, double t0_inj,
                       double ts_inj,
                       double E0_refresh, double k_refresh, double umin_refresh,
                       double thetaCore, double theta0,
                       double thetaCoreGlobal);

double shockVel(double u);
double E_inj(double te, double L0, double q, double ts);
double L_inj(double te, double L0, double q, double ts);
double envDensityPar(double R, struct shockParams *pars);
double envMassPar(double R, struct shockParams *pars);
double envRadiusPar(double M, struct shockParams *pars);
double envDensity(double R, int envType, double rho0, double R0, double k,
                  double rho1);
double envMass(double R, int envType, double rho0, double R0, double k,
               double rho1);
double envRadius(double M, int envType, double rho0, double R0, double k,
                 double rho1);
void shockInitDecel(double t0, double *R0, double *u0, void *argv); 
void shockInitFind(double t0, double *R0, double *u0, double tRes, void *argv);
void shockInitFindISM(double t0, double *R0, double *u0, double tRes,
                      void *argv);
void Rudot2D(double t, double *x, void *argv, double *xdot);
void RuThdot3D(double t, double *x, void *argv, double *xdot);
void shockEvolveRK4(double *t, double *R, double *u, int N, double R0, 
                    double u0, void *args);
void shockEvolveSpreadRK4(double *t, double *R, double *u, double *th, int N, 
                            double R0, double u0, double th0, void *args);

#endif
