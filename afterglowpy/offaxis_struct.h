#ifndef AFTERGLOWPY_STRUCT
#define AFTERGLOWPY_STRUCT

// offaxis.h

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "integrate.h"

#define ERR_CHK_VOID(pars) if(pars->error){ return;}
#define ERR_CHK_INT(pars) if(pars->error){ return 0;}
#define ERR_CHK_DBL(pars) if(pars->error){ return 0.0;}
#define MSG_LEN 4096
#define DUMP_MSG_LEN_MAX 16384  //overkill: 200 lines * 80c per line = 16000

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// some physical and mathematical constants
#define PI          3.14159265358979323846
#define v_light     2.99792458e10   // speed of light in cm / s
#define invv_light  3.335640952e-11 // inverse speed of light s / cm
#define m_e         9.1093897e-28   // electron mass in g
#define m_p         1.6726231e-24   // proton mass in g
#define invm_e      1.097768383e27  // inverse electron mass in 1/g
#define invm_p      5.978633202e23  // inverse proton mass in 1/g
#define h_planck    6.6260755e-27   // Planck's constant in erg / s
#define h_bar       1.05457266e-27  // Planck's constant / 2 PI in erg /s
#define k_B         1.380658e-16    // Boltzmann's constant in erg / K
#define e_e         4.803e-10       // electron charge in Gaussian cgs units
#define sigma_T     6.65e-25        // Thomson cross section free electron cm^2
#define cgs2mJy     1e26            // quantity in cgs to mJy
#define mJy2cgs     1e-26           // quantity in mJy to cgs
#define deg2rad     0.017453292     // quantity in degrees to radians
#define rad2deg     57.29577951     // quantity in radians to degrees
#define sec2day     0.000011574     // quantity in seconds to days
#define day2sec     86400           // quantity in days to seconds
#define parsec      3.0857e18       // quantity in parsec to cm
#define Hz2eV       4.13566553853599E-15
#define eV2Hz       2.417991e+14

#define _cone -2
#define _tophat -1
#define _Gaussian 0
#define _powerlaw_core 1 //has a core as well
#define _Gaussian_core 2 // has a core as well
#define _spherical 3
#define _powerlaw 4
#define _exponential 5
#define _twocomponent 6
#define _exponential2 7

#define IC_COOLING_FLAG 1
#define EPS_E_BAR_FLAG  2
#define SSA_SMOOTH_FLAG 4
#define SSA_SHARP_FLAG  8

enum{INT_TRAP_FIXED, INT_TRAP_ADAPT, INT_SIMP_FIXED, INT_SIMP_ADAPT,
     INT_ROMB_ADAPT, INT_TRAP_NL, INT_HYBRID, INT_CADRE,
     INT_GK49_ADAPT, INT_GK715_ADAPT, INT_GK1021_ADAPT,
     INT_UNDEFINED};

enum{GAMMA_INF, GAMMA_FLAT, GAMMA_EVENMASS, GAMMA_STRUCT};

struct fluxParams
{
    double theta;
    double phi;
    double cp;
    double ct;
    double st;
    double cto;
    double sto;
    
    double theta_obs;
    double t_obs;
    double nu_obs;
    double d_L;

    double E_iso;
    double n_0;
    double g_init;

    double p;
    double epsilon_E;
    double epsilon_B;
    double ksi_N;

    double theta_h;
    double E_iso_core;
    double theta_core;
    double theta_wing;
    double b;
    double E_tot;
    double g_core;
    double E_core_global;
    double theta_core_global;

    int envType;
    double As;
    double Rwind;

    double L0;
    double q;
    double ts;
    
    double current_theta_cone_hi;
    double current_theta_cone_low;
    double theta_obs_cur;
    int tRes;
    int latRes;
    int spread;
    int counterjet;

    int int_type;
    double rtol_struct;
    double rtol_theta;
    double rtol_phi;
    int nmax_theta;
    int nmax_phi;

    double atol_theta;

    double Rt0;
    double Rt1;
    double ta;
    double tb;

    double C_BMsqrd;
    double C_STsqrd;

    double t_NR;

    double *t_table;
    double *R_table;
    double *u_table;
    double *th_table;
    double *mu_table;
    int table_entries;

    double *t_table_inner;
    double *R_table_inner;
    double *u_table_inner;
    double *th_table_inner;
    double *mu_table_inner;
    int table_entries_inner;

    int spec_type;
    int gamma_type;

    double (*f_E)(double, void *);

    double *mask;
    int nmask;

    long nevals;

    int error;
    char *error_msg;
};


double dmin(const double a, const double b);


double f_E_tophat(double theta, void *params);
double f_E_Gaussian(double theta, void *params);
double f_E_powerlaw(double theta, void *params);
double f_E_twocomponent(double theta, void *params);
double f_E_exponential(double theta, void *params);
double f_Etot_tophat(void *params);
double f_Etot_Gaussian(void *params);
double f_Etot_powerlaw(void *params);

void make_R_table(struct fluxParams *pars);
void make_mu_table(struct fluxParams *pars);
double check_t_e(double t_e, double mu, double t_obs, double *mu_table, int N);
int searchSorted(double x, double *arr, int N);
double interpolateLin(int a, int b, double x, double *X, double *Y, int N);
double interpolateLog(int a, int b, double x, double *X, double *Y, int N);
double find_jet_edge(double phi, double cto, double sto, double theta0,
                     double *a_mu, double *a_thj, int N);
double costheta_integrand(double a_theta, void* params); // inner integral
double phi_integrand(double a_phi, void* params); // outer integral
double emissivity(double nu, double R, double mu, double te,
                    double u, double us, double n0, double p, double epse,
                    double epsB, double ksiN, int specType); //emissivity of
                                                             // a zone.
double flux(struct fluxParams *pars, double atol); // determine flux for a given t_obs

double flux_cone(double t_obs, double nu_obs, double E_iso, double theta_h,
                    double theta_cone_low, double theta_cone_hi,
                    double atol, struct fluxParams *pars);
double intensity(double theta, double phi, double tobs, double nuobs,
                double theta_obs, double theta_cone_hi, double theta_cone_low,
                struct fluxParams *pars);
void shockVals(double theta, double phi, double tobs,
                 double *t, double *R, double *u, double *thj,
                 double theta_obs, double theta_cone_hi, double theta_cone_low,
                 struct fluxParams *pars);
void intensity_cone(double *theta, double *phi, double *t, double *nu, 
                        double *I, int N, double E_iso_core, 
                        double theta_h_core, double theta_h_wing, 
                        struct fluxParams *pars);
void intensity_struct(double *theta, double *phi, double *t, double *nu, 
                        double *I, int N,
                        double E_iso_core, 
                        double theta_h_core, double theta_h_wing,
                        int res_cones, double (*f_E)(double,void *),
                        struct fluxParams *pars);
void intensity_structCore(double *theta, double *phi, double *t, double *nu, 
                        double *I, int N,
                        double E_iso_core, 
                        double theta_h_core, double theta_h_wing,
                        int res_cones, double (*f_E)(double,void *),
                        struct fluxParams *pars);
void shockVals_cone(double *theta, double *phi, double *tobs, 
                   double *t, double *R, double *u, double *thj, int N,
                   double E_iso_core, double theta_h_core, double theta_h_wing, 
                   struct fluxParams *pars);
void shockVals_struct(double *theta, double *phi, double *tobs,
                        double *t, double *R, double *u, double *thj, int N,
                        double E_iso_core, 
                        double theta_h_core, double theta_h_wing,
                        int res_cones, double (*f_E)(double,void *),
                        struct fluxParams *pars);
void shockVals_structCore(double *theta, double *phi, double *tobs,
                        double *t, double *R, double *u, double *thj, int N,
                        double E_iso_core, 
                        double theta_h_core, double theta_h_wing,
                        int res_cones, double (*f_E)(double,void *),
                        struct fluxParams *pars);
void lc_tophat(double *t, double *nu, double *F, int Nt,
                double E_iso, double theta_h, struct fluxParams *pars);
void lc_cone(double *t, double *nu, double *F, int Nt, double E_iso,
                double theta_h, double theta_wing, struct fluxParams *pars);
void lc_powerlawCore(double *t, double *nu, double *F, int Nt,
                    double E_iso_core, double theta_h_core, 
                    double theta_h_wing, double beta,
                    double *theta_c_arr, double *E_iso_arr,
                    int res_cones, struct fluxParams *pars);
void lc_powerlaw(double *t, double *nu, double *F, int Nt,
                    double E_iso_core, double theta_h_core, 
                    double theta_h_wing,
                    double *theta_c_arr, double *E_iso_arr,
                    int res_cones, struct fluxParams *pars);
void lc_Gaussian(double *t, double *nu, double *F, int Nt,
                        double E_iso_core, 
                        double theta_h_core, double theta_h_wing,
                        double *theta_c_arr, double *E_iso_arr,
                        int res_cones, struct fluxParams *pars);
void lc_GaussianCore(double *t, double *nu, double *F, int Nt,
                        double E_iso_core,
                        double theta_h_core, double theta_h_wing,
                        double *theta_c_arr, double *E_iso_arr,
                        int res_cones, struct fluxParams *pars);

void calc_flux_density(int jet_type, int spec_type, 
                            double *t, double *nu, double *Fnu, int N,
                            struct fluxParams *fp);
void calc_intensity(int jet_type, int spec_type, double *theta, double *phi,
                            double *t, double *nu, double *Inu, int N,
                            struct fluxParams *fp);
void calc_shockVals(int jet_type, double *theta, double *phi, double *tobs,
                            double *t, double *R, double *u, double *thj, int N,
                            struct fluxParams *fp);

void setup_fluxParams(struct fluxParams *pars,
                    double d_L,
                    double theta_obs,
                    double E_iso_core, double theta_core, double theta_wing,
                    double b, double L0, double q, double ts, 
                    double n_0,
                    double p,
                    double epsilon_E,
                    double epsilon_B, 
                    double ksi_N,
                    double g0,
                    double E_core_global,
                    double theta_core_global,
                    double ta, double tb,
                    int tRes, int latRes, int int_type,
                    double rtol_struct, double rtol_phi, double rtol_theta,
                    int nmax_phi, int nmax_theta,
                    int spec_type,
                    double *mask, int nmask,
                    int spread, int counterjet, int gamma_type);
void set_jet_params(struct fluxParams *pars, double E_iso, double theta_h);
void set_obs_params(struct fluxParams *pars, double t_obs, double nu_obs,
                        double theta_obs_cur, double current_theta_cone_hi, 
                        double current_theta_cone_low);
int check_error(void *params);
void set_error(struct fluxParams *pars, char msg[]);
void free_fluxParams(struct fluxParams *pars);

#endif
