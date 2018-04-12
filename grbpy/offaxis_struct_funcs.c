#include "offaxis_struct.h"

double dmin(const double a, const double b)
{
    if(a <= b)
        return a;
    else
        return b;
}
 
/////////////////////////////////////////////////////////////////////////

double f_E_tophat(double theta, void *params)
{
    struct fluxParams *pars = (struct fluxParams *)params;
    if(theta <= pars->theta_core)
        return pars->E_iso_core;
    return 0.0;
}

double f_E_Gaussian(double theta, void *params)
{
    struct fluxParams *pars = (struct fluxParams *)params;
    if(theta <= pars->theta_wing)
    {
        double x = theta / pars->theta_core;
        return pars->E_iso_core * exp(-0.5*x*x);
    }
    return 0.0;
}

double f_E_powerlaw(double theta, void *params)
{
    struct fluxParams *pars = (struct fluxParams *)params;
    if(theta <= pars->theta_wing)
    {
        double x = theta / pars->theta_core;
        return pars->E_iso_core / (1 + x*x);
    }
    return 0.0;
}

double f_Etot_tophat(void *params)
{
    struct fluxParams *pars = (struct fluxParams *)params;
    double E0 = pars->E_iso_core;
    double thetaC = pars->theta_core;
    return E0*(1-cos(thetaC));
}

double f_Etot_Gaussian(void *params)
{
    struct fluxParams *pars = (struct fluxParams *)params;
    double E0 = pars->E_iso_core;
    double thC = pars->theta_core;
    double thW = pars->theta_wing;
    double a = 0.5*thW*thW/(thC*thC);
    double expa = exp(-a);
    double I0 = 1-expa;
    double I1 = 1-(1 + a)*expa;
    double I2 = 2-(2 + 2*a + a*a)*expa;
    double I3 = 6-(6 + 6*a + 3*a*a + a*a*a)*expa;
    return E0*thC*thC*(I0 - thC*thC*I1/3.0 + thC*thC*thC*thC*I2/30.0
                        - thC*thC*thC*thC*thC*thC*I3/630);
}

double f_Etot_powerlaw(void *params)
{
    struct fluxParams *pars = (struct fluxParams *)params;
    double E0 = pars->E_iso_core;
    double thC = pars->theta_core;
    double thW = pars->theta_wing;
    double a = thW*thW/(thC*thC);
    double atana = atan(a);
    double la = log(a*a+1);
    double I0 = atana;
    double I1 = 0.5*la;
    double I2 = a - atana;
    double I3 = 0.5*a*a - 0.5*la;
    return 0.5*E0*thC*thC*(I0 - thC*thC*I1/3.0 + thC*thC*thC*thC*I2/30.0
                        - thC*thC*thC*thC*thC*thC*I3/630);
}

///////////////////////////////////////////////////////////////////////////////

double get_lfacbetashocksqrd(double a_t_e, double C_BMsqrd, double C_STsqrd)
{
    return C_BMsqrd / (a_t_e*a_t_e*a_t_e) + C_STsqrd * pow(a_t_e, -6.0/5.0);
}

///////////////////////////////////////////////////////////////////////////////

double get_lfacbetasqrd(double a_t_e, double C_BMsqrd, double C_STsqrd)
{
    return 0.5 * C_BMsqrd / (a_t_e*a_t_e*a_t_e) + 9.0 / 16.0 * C_STsqrd * 
        pow(a_t_e, -6.0/5.0);
}

///////////////////////////////////////////////////////////////////////////////

double check_t_e(double t_e, double mu, double t_obs, double *mu_table, int N)
{
    if(mu > mu_table[N-1])
    {
        printf("mu >> 1? this should not have happened\n");
        printf("   t_obs=%.6lg t_e=%.6lg mu=%.6lg mu_table[-1]=%.6lg\n",
                t_obs, t_e, mu, mu_table[N-1]);
        abort();
    }

    if(mu_table[0] >= mu) // happens only if t_e very small
    {
        printf("very small mu: mu=%.3lg, mu[0]=%.3lg\n", mu, mu_table[0]);
        return t_obs / (1.0 - mu); // so return small t_e limit
    }

    return t_e;
}

int searchSorted(double x, double *arr, int N)
{
    if(x <= arr[0])
        return 0;
    else if(x >= arr[N-1])
        return N-2;

  unsigned int i = ((unsigned int) N) >> 1;
  unsigned int a = 0;
  unsigned int b = N-1;
  
  while (b-a > 1u)
  {
    i = (b+a) >> 1;
    if (arr[i] > x)
        b = i;
    else
        a = i;
  }

  return (int)a;
}

double interpolateLin(int a, int b, double x, double *X, double *Y, int N)
{
    double xa = X[a];
    double xb = X[b];
    double ya = Y[a];
    double yb = Y[b];

    return ya + (yb-ya) * (x-xa)/(xb-xa);
}

double interpolateLog(int a, int b, double x, double *X, double *Y, int N)
{
    double xa = X[a];
    double xb = X[b];
    double ya = Y[a];
    double yb = Y[b];

    return ya * pow(yb/ya, log(x/xa)/log(xb/xa));
}
///////////////////////////////////////////////////////////////////////////////

double Rintegrand(double a_t_e, void* params)
{
    double C_BMsqrd = ((double *) params)[0];
    double C_STsqrd = ((double *) params)[1];
    double us2 = get_lfacbetashocksqrd(a_t_e, C_BMsqrd, C_STsqrd);
    if(!isfinite(us2))
        return v_light;
    return v_light * sqrt(us2 / (1.0 + us2));
}

///////////////////////////////////////////////////////////////////////////////


void make_mu_table(struct fluxParams *pars)
{
    double t_obs = pars->t_obs;
    double *t_table = pars->t_table;
    double *R_table = pars->R_table;
    double *mu_table = pars->mu_table;
    int table_entries = pars->table_entries;

    int i;
    for (i = 0; i< table_entries; i++)
    {
        mu_table[i] = (t_table[i] - t_obs) / R_table[i] * v_light;
    }
}

void make_R_table(struct fluxParams *pars)
{
    // min/max observer times
    double ta = pars->ta;
    double tb = pars->tb;

    double C_BMsqrd = pars->C_BMsqrd;
    double C_STsqrd = pars->C_STsqrd;
    double t_NR = pars->t_NR;

    double Rt0, Rt1;
    //Rt0 = 1.0e-2 * day2sec;
    //Rt1 = 1.0e7 * day2sec; 

    //at fixed t_obs, earliest emission is *always* from mu=-1
    // so t_obs ~ t_e
    Rt0 = 0.1*ta;

    //at fixed t_obs, latest emission is *always* from mu=+1
    // so t_obs ~ t-R/c
    if(tb > 0.1*t_NR)  // at late times R<<c so t_obs ~ t_e
        Rt1 = 10*tb;
    else // at early times t_obs ~ t*(gamma_sh^-2)/8 ~ CBM^-2 * t^4 / 8
        Rt1 = 10*pow(8*tb*C_BMsqrd, 0.25);
    
    //printf("        Rt0: %.1le Rt1: %.1le\n", Rt0, Rt1);

    int tRes = pars->tRes;
    int table_entries = (int)(tRes * log10(Rt1/Rt0));
    
    pars->table_entries = table_entries;
    pars->Rt0 = Rt0;
    pars->Rt1 = Rt1;

    pars->t_table = (double *)realloc(pars->t_table, 
                                        sizeof(double) * table_entries);
    pars->R_table = (double *)realloc(pars->R_table, 
                                        sizeof(double) * table_entries);
    pars->mu_table = (double *)realloc(pars->mu_table, 
                                        sizeof(double) * table_entries);
    double *t_table = pars->t_table;
    double *R_table = pars->R_table;

    double DR, R;
    double t;
    double tp = 0.0; // time for previous table entry
    double Rp = 0.0; // R value of previous table entry
    int i;

    // prepare integration function
    double Rpar[2] = {C_BMsqrd, C_STsqrd};
#ifdef USEGSL
    gsl_function F;
    F.function = &Rintegrand;
    F.params = Rpar;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    double error;
#endif

    // set up R table. Each entry is equal to previous plus additional distance
    double fac0 = pow(Rt1/Rt0, 1.0/(table_entries-1.0));
    double fac = 1.0;
    for (i=0; i < table_entries-1; i++)
    {
        //t = Rt0 * pow( Rt1 / Rt0, (double) i / (double) (table_entries - 1.0));
        t = Rt0 * fac;
        fac *= fac0;
#ifdef USEGSL
        gsl_integration_qag (&F, tp, t, 0, 1.0e-6, 1000, 1, w, &DR, &error);
#else
        DR = romb(&Rintegrand, tp, t, 1000, 0, R_ACC, Rpar);
#endif
        R = Rp + DR;
        t_table[i] = t; R_table[i] = R;
        Rp = R; tp = t;
    }
    t = Rt1;
#ifdef USEGSL
    gsl_integration_qag (&F, tp, t, 0, 1.0e-6, 1000, 1, w, &DR, &error);
#else
    DR = romb(&Rintegrand, tp, t, 1000, 0, R_ACC, Rpar);
#endif
    R = Rp + DR;
    t_table[table_entries-1] = t;
    R_table[table_entries-1] = R;
    Rp = R; tp = t;

    if(R_table[0] != R_table[0])
    {
        printf("Rintegration Error: R[0]=%.3e  (fac0=%.3e)\n", 
                                    R_table[0], fac0);
        printf("    t=%.3e Rdot=%.3e t=%.3e Rdot=%.3e\n",
                    0.0,Rintegrand(0.0,Rpar), t_table[0], 
                    Rintegrand(t_table[0],Rpar));
    }

    // free memory for integration routine
#ifdef USEGSL
    gsl_integration_workspace_free(w);
#endif
}

///////////////////////////////////////////////////////////////////////////////

double emissivity(double nu, double R, double sinTheta, double mu, double te,
                    double u, double us, double n0, double p, double epse,
                    double epsB, double ksiN, int specType)
{
    if(us < 1.0e-5)
    {
        //shock is ~ at sound speed of warm ISM. Won't shock, approach invalid.
        return 0.0;
    }
    if(sinTheta == 0.0 || R == 0.0)
        return 0.0;

    // set remaining fluid quantities
    double g = sqrt(1+u*u);
    double beta = u/g;
    double betas = us / sqrt(1+us*us);
    double nprime = 4.0 * n0 * g; // comoving number density
    double e_th = u*u/(g+1) * nprime * m_p * v_light * v_light;
    double B = sqrt(epsB * 8.0 * PI * e_th);
    double a = (1.0 - mu * beta); // beaming factor
    double ashock = (1.0 - mu * betas); // shock velocity beaming factor
    double DR = R / (12.0 * g*g * ashock);
    if (DR < 0.0) DR *= -1.0; // DR is function of the absolute value of mu


    // set local emissivity 
    double nuprime = nu * g * a; // comoving observer frequency
    double g_m = (2.0 - p) / (1.0 - p) * epse * e_th / (
                            ksiN * nprime * m_e * v_light * v_light);
    double g_c = 6 * PI * m_e * g * v_light / (sigma_T * B * B * te);

    //Inverse Compton adjustment of lfac_c
    if(specType == 1)
    {
        double gr = g_c / g_m;
        double y = beta * epse/epsB;
        double X = 1.0;

        if(gr <= 1.0 || gr*gr-gr-y <= 0.0)
        {
            //Fast Cooling
            X = 0.5*(1 + sqrt(1+4*y));
        }
        else
        {
            //Slow Cooling
            double b = y * pow(gr, 2-p);
            double Xa = 1 + b;
            double Xb = pow(b, 1.0/(4-p)) + 1.0/(4-p);
            double s = b*b / (b*b + 1);
            X = Xa * pow(Xb/Xa, s);
            int i;
            for(i=0; i<5; i++)
            {
                double po = pow(X, p-2);
                double f = X*X - X - b*po;
                double df = 2*X - 1 - (p-2)*b*po/X;
                double dX = -f/df;
                X += dX;
                if(fabs(dX) < 1.0e-4*X)
                break;
            }
        }

        g_c /= X;
    }

    double nu_m = 3.0 * g_m * g_m * e_e * B / (4.0 * PI * m_e * v_light);
    double nu_c = 3.0 * g_c * g_c * e_e * B / (4.0 * PI * m_e * v_light);
    double em = 0.5*(p - 1.0)*sqrt(3.0) * e_e*e_e*e_e * ksiN * nprime * B
                    / (m_e*v_light*v_light);
  
    double freq = 0.0; // frequency dependent part of emissivity


    // set frequency dependence
    if (nu_c > nu_m)
    {
        if (nuprime < nu_m) 
            freq = pow(nuprime / nu_m, 1.0 / 3.0 );
        else if (nuprime < nu_c)
            freq = pow(nuprime / nu_m, 0.5 * (1.0 - p));
        else
            freq = pow(nu_c / nu_m, 0.5 * (1.0 - p))
                    * pow(nuprime / nu_c, -0.5*p);
    }
    else
    {
        if (nuprime < nu_c)
            freq = pow(nuprime / nu_c, 1.0/3.0);
        else if (nuprime < nu_m)
            freq = sqrt(nu_c / nuprime);
        else
            freq = sqrt(nu_c/nu_m) * pow(nuprime / nu_m, -0.5 * p);
    }


    if(em != em || em < 0.0)
        printf("bad em at:%.3le te=%.3le sinTheta=%.3lf mu=%.3lf\n",
                em, te, sinTheta, mu);
    if(freq != freq || freq < 0.0)
        printf("bad freq at:%.3le te=%.3le sinTheta=%.3lf mu=%.3lf\n",
                freq, te, sinTheta, mu);

    return R * R * sinTheta * DR * em * freq / (g*g * a*a);
}

double theta_integrand(double a_theta, void* params) // inner integral
{
    struct fluxParams *pars = (struct fluxParams *) params;
    
    //double cp = cos(pars->phi); 
    //double cto = cos(pars->theta_obs_cur);
    //double sto = sin(pars->theta_obs_cur);
    double ast = sin(a_theta);
    double act = cos(a_theta);
    double mu = ast * (pars->cp) * (pars->sto) + act * (pars->cto);

    int ia = searchSorted(mu, pars->mu_table, pars->table_entries);
    int ib = ia+1;
    double t_e = interpolateLin(ia, ib, mu, pars->mu_table, pars->t_table, 
                            pars->table_entries);
    t_e = check_t_e(t_e, mu, pars->t_obs, pars->mu_table, pars->table_entries);
    
    double R = interpolateLog(ia, ib, t_e, pars->t_table, pars->R_table, 
                            pars->table_entries);

    //printf("%e, %e, %e # tobs, R, t_e\n", t_obs, t_e, R);
    double us2 = get_lfacbetashocksqrd(t_e, pars->C_BMsqrd, 
                                                    pars->C_STsqrd);
    double u2 = get_lfacbetasqrd(t_e, pars->C_BMsqrd, pars->C_STsqrd);
    
    double dFnu =  emissivity(pars->nu_obs, R, ast, mu, t_e, sqrt(u2), 
                                sqrt(us2), pars->n_0, pars->p, pars->epsilon_E,
                                pars->epsilon_B, pars->ksi_N, pars->spec_type);

    int i;
    double fac = 1.0;
    for(i=0; i<pars->nmask; i++)
    {
        double *m = &((pars->mask)[9*i]);
        if(m[0]<t_e && t_e<m[1] && m[2]<R && R<m[3] && m[4]<a_theta
                && a_theta<m[5] && m[6]<pars->phi && pars->phi<m[7])
            fac = m[8];
    }

    return fac * dFnu;
}

///////////////////////////////////////////////////////////////////////////////

double phi_integrand(double a_phi, void* params) // outer integral
{
  double result;

  struct fluxParams *pars = (struct fluxParams *) params;
  
  // set up integration routine
#ifdef USEGSL
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  gsl_function F; F.function = &theta_integrand; F.params = params;
  double error;
#endif

  pars->phi = a_phi;
  pars->cp = cos(a_phi);
  
  // implement sideways spreading approximation until spherical symmetry reached
  double theta_1 = pars->current_theta_cone_hi;
  double theta_0 = pars->current_theta_cone_low;
  double Dtheta = theta_1 - theta_0;
  if (pars->t_obs > pars->t_NR)
  {
    theta_1 = dmin(0.5 * PI, 
                    pars->theta_h + 0.1 * log( pars->t_obs / pars->t_NR));
  }
 
  //printf("# theta integration domain: %e - %e\n", theta_1 - Dtheta, theta_1); fflush(stdout);
 
  // For a given phi, integrate over theta
#ifdef USEGSL
  gsl_integration_qags(&F, theta_1 - Dtheta, theta_1, 0, 1.0e-4, 1000, w, 
                        &result, &error);
  // free integration routine memory
  gsl_integration_workspace_free(w);
#else
  result = romb(&theta_integrand, theta_1-Dtheta, theta_1, 1000, 
                    pars->theta_atol, THETA_ACC, params);
#endif
  if(result != result || result < 0.0)
      printf("bad result at:%.3le t_obs=%.3le theta_lo=%.3lf theta_hi=%.3lf phi=%.3lf\n",
              result, pars->t_obs, theta_0, theta_1, pars->phi);
  
  //return result
  return result;
}

double phi_integrand_vec(double phi, void *params)
{
    struct fluxParams *pars = (struct fluxParams *) params;
    
    double cp = cos(phi); 
    double mu = cp * (pars->st) * (pars->sto) + (pars->ct) * (pars->cto);

    int ia = searchSorted(mu, pars->mu_table, pars->table_entries);
    int ib = ia+1;
    double t_e = interpolateLin(ia, ib, mu, pars->mu_table, pars->t_table, 
                                pars->table_entries);
    t_e = check_t_e(t_e, mu, pars->t_obs, pars->mu_table, pars->table_entries);
    
    double R = interpolateLog(ia, ib, t_e, pars->t_table, pars->R_table, 
                            pars->table_entries);

    //printf("%e, %e, %e # tobs, R, t_e\n", t_obs, t_e, R);
    double us2 = get_lfacbetashocksqrd(t_e, pars->C_BMsqrd, 
                                                    pars->C_STsqrd);
    double u2 = get_lfacbetasqrd(t_e, pars->C_BMsqrd, pars->C_STsqrd);
    
    double dFnu =  emissivity(pars->nu_obs, R, pars->st, mu, t_e, sqrt(u2), 
                                sqrt(us2), pars->n_0, pars->p, pars->epsilon_E,
                                pars->epsilon_B, pars->ksi_N, pars->spec_type);

    return dFnu;
}

void theta_integrand_vec(double theta, double *Fnu, double *t, double *nu,
                            int Nt, void *params)
{
    struct fluxParams *pars = (struct fluxParams *)params;

    double E = pars->f_E(theta, params);
    set_jet_params(pars, E, theta);

    pars->theta = theta;
    pars->ct = cos(theta);
    pars->st = sin(theta);

    int i;
    for(i=0; i<Nt; i++)
    {
        double theta_obs = pars->theta_obs;
        set_obs_params(pars, t[i], nu[i], theta_obs, theta, theta);
        make_mu_table(pars);
        double F1 = 2.0 * romb(&phi_integrand_vec, 0.0, PI, 1000, 0, PHI_ACC,
                                                params);

        //Counter-jet
        theta_obs = PI - pars->theta_obs;
        set_obs_params(pars, t[i], nu[i], theta_obs, theta, theta);
        double F2 = 2.0 * romb(&phi_integrand_vec, 0.0, PI, 1000, 0, PHI_ACC,
                                                params);
        Fnu[i] = F1 + F2;
    }
}

///////////////////////////////////////////////////////////////////////////////

double flux(struct fluxParams *pars, double atol) // determine flux for a given t_obs
{
  double result;
  double phi_0 = 0.0;
  double phi_1 = PI;
  
  // set up integration routines for integration over phi
#ifdef USEGSL
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F;
    F.function = &phi_integrand;
    F.params = pars;
    double error;
#endif

  // at this stage t_obs is known, so mu_table can be made
  make_mu_table(pars); 

  double d_L = pars->d_L;

  double Fcoeff = cgs2mJy / (4*PI * d_L*d_L);
  
  //printf("about to integrate phi between %e and %e\n", phi_0, phi_1); fflush(stdout);
#ifdef USEGSL
  gsl_integration_qags (&F, phi_0, phi_1, 0, 1.0e-3, 1000, w, 
                            &result, &error); 
  // free memory
  gsl_integration_workspace_free(w);
#else
  //pars->theta_atol = 0.0;
  //double I0 = phi_integrand(0.0, pars);
  //pars->theta_atol = 1.0e-6 * I0;
  result = 2 * Fcoeff * romb(&phi_integrand, phi_0, phi_1, 1000, 
                                            atol/(2*Fcoeff), PHI_ACC, pars);
#endif

  //return result

  return result;
}

void lc_cone(double *t, double *nu, double *F, int Nt, double E_iso,
                double theta_core, double theta_wing, struct fluxParams *pars)
{
    int i;

    set_jet_params(pars, E_iso, 0.5*(theta_core+theta_wing));

    for(i=0; i<Nt; i++)
        F[i] = flux_cone(t[i], nu[i], -1, -1, theta_core, theta_wing, 0.0,
                            pars);
}

void lc_tophat(double *t, double *nu, double *F, int Nt,
                double E_iso, double theta_h, struct fluxParams *pars)
{
    int i;

    set_jet_params(pars, E_iso, theta_h);

    for(i=0; i<Nt; i++)
        F[i] = flux_cone(t[i], nu[i], -1, -1, 0.0, theta_h, 0.0, pars);
}

void lc_powerlawCore(double *t, double *nu, double *F, int Nt, 
                    double E_iso_core, double theta_h_core, 
                    double theta_h_wing, double beta,
                    double *theta_c_arr, double *E_iso_arr,
                    int res_cones, struct fluxParams *pars)
{
    //Flux from a powerlaw jet, Eiso ~ theta^beta

    //Core
    lc_tophat(t, nu, F, Nt, E_iso_core, theta_h_core, pars);

    double Dtheta, theta_cone_hi, theta_cone_low, theta_h, theta_c, E_iso;

    Dtheta = (theta_h_wing - theta_h_core) / res_cones;

    int i,j;
    for(i=0; i<res_cones; i++)
    {
        theta_c = theta_h_core + (i+0.5)*Dtheta;
        E_iso = E_iso_core * pow(theta_c/theta_h_core, beta);

        theta_cone_hi = theta_h_core + (i+1)*Dtheta;
        theta_cone_low = theta_h_core + i*Dtheta;
        theta_h = theta_cone_hi;
    
        if(theta_c_arr != NULL)
            theta_c_arr[i] = theta_c;
        if(E_iso_arr != NULL)
            E_iso_arr[i] = E_iso;

        set_jet_params(pars, E_iso, theta_h);

        for(j=0; j<Nt; j++)
            F[j] += flux_cone(t[j], nu[j], -1, -1, theta_cone_low, 
                                theta_cone_hi, F[j]*pars->flux_rtol/res_cones,
                                pars);
    }
}

void lc_powerlaw(double *t, double *nu, double *F, int Nt,
                        double E_iso_core, 
                        double theta_h_core, double theta_h_wing,
                        double *theta_c_arr, double *E_iso_arr,
                        int res_cones, struct fluxParams *pars)
{
    //Flux from a smooth powerlaw
    
    int i,j;
    //No Core
    for(j=0; j<Nt; j++)
        F[j] = 0.0;

    double Dtheta, theta_cone_hi, theta_cone_low, theta_h, theta_c, E_iso;

    //if(theta_h_wing > 10*theta_h_core)
    //        theta_h_wing = 10*theta_h_core;

    Dtheta = theta_h_wing / res_cones;

    for(i=0; i<res_cones; i++)
    {
        theta_c = (i+0.5) * Dtheta;
        E_iso = E_iso_core
                    / (1.0 + theta_c*theta_c/(theta_h_core*theta_h_core));

        theta_cone_hi = (i+1) * Dtheta;
        theta_cone_low = i * Dtheta;
        theta_h = theta_cone_hi;

        if(theta_c_arr != NULL)
            theta_c_arr[i] = theta_c;
        if(E_iso_arr != NULL)
            E_iso_arr[i] = E_iso;


        set_jet_params(pars, E_iso, theta_h);

        for(j=0; j<Nt; j++)
            F[j] += flux_cone(t[j], nu[j], -1, -1, theta_cone_low,
                                theta_cone_hi, F[j]*pars->flux_rtol/res_cones,
                                pars);
    }
}

void lc_Gaussian(double *t, double *nu, double *F, int Nt,
                        double E_iso_core, 
                        double theta_h_core, double theta_h_wing,
                        double *theta_c_arr, double *E_iso_arr,
                        int res_cones, struct fluxParams *pars)
{
    //Flux from a Gaussian jet.
    
    int i,j;
    //No Core
    for(j=0; j<Nt; j++)
        F[j] = 0.0;

    pars->E_tot = f_Etot_Gaussian(pars);

    double Dtheta, theta_cone_hi, theta_cone_low, theta_h, theta_c, E_iso;

    //if(theta_h_wing > 10*theta_h_core)
    //    theta_h_wing = 10*theta_h_core;

    Dtheta = theta_h_wing / res_cones;

    for(i=0; i<res_cones; i++)
    {
        //printf("Gaussian: thetaC = %.6f\n", theta_c);
        theta_c = (i+0.5) * Dtheta;
        E_iso = E_iso_core
                    * exp(-0.5 * theta_c*theta_c/(theta_h_core*theta_h_core));

        theta_cone_hi = (i+1) * Dtheta;
        theta_cone_low = i * Dtheta;
        theta_h = theta_cone_hi;

        if(theta_c_arr != NULL)
            theta_c_arr[i] = theta_c;
        if(E_iso_arr != NULL)
            E_iso_arr[i] = E_iso;


        set_jet_params(pars, E_iso, theta_h);

        for(j=0; j<Nt; j++)
            F[j] += flux_cone(t[j], nu[j], -1, -1, theta_cone_low,
                                theta_cone_hi, F[j]*pars->flux_rtol/res_cones,
                                pars);
    }
}

void lc_GaussianCore(double *t, double *nu, double *F, int Nt,
                        double E_iso_core,
                        double theta_h_core, double theta_h_wing,
                        double *theta_c_arr, double *E_iso_arr,
                        int res_cones, struct fluxParams *pars)
{
    //Flux from a Gaussian jet, with a core.

    //Core
    lc_tophat(t, nu, F, Nt, E_iso_core, theta_h_core, pars);

    double Dtheta, theta_cone_hi, theta_cone_low, theta_h, theta_c, E_iso;

    //if(theta_h_wing > 10*theta_h_core)
    //    theta_h_wing = 10*theta_h_core;

    Dtheta = (theta_h_wing - theta_h_core) / res_cones;

    int i, j;
    for(i=0; i<res_cones; i++)
    {
        theta_c = theta_h_core + (i+0.5) * Dtheta;
        E_iso = E_iso_core * exp(-0.5
                * (theta_c*theta_c/(theta_h_core*theta_h_core) - 1.0));

        theta_cone_hi = theta_h_core + (i+1) * Dtheta;
        theta_cone_low = theta_h_core + i * Dtheta;
        theta_h = theta_cone_hi;

        if(theta_c_arr != NULL)
            theta_c_arr[i] = theta_c;
        if(E_iso_arr != NULL)
            E_iso_arr[i] = E_iso;

        set_jet_params(pars, E_iso, theta_h);

        for(j=0; j<Nt; j++)
            F[j] += flux_cone(t[j], nu[j], -1, -1, theta_cone_low,
                                theta_cone_hi, F[j]*pars->flux_rtol/res_cones,
                                pars);
    }
}

void lc_vec(double *t, double *nu, double *Fnu, int Nt, double E_iso_core,
            double theta_core, double theta_wing, int Ntheta, 
            double (*f_E)(double, void *), double (*f_Etot)(void *),
            struct fluxParams *pars)
{
    pars->E_iso_core = E_iso_core;
    pars->theta_core = theta_core;
    pars->theta_wing = theta_wing;
    pars->f_E = f_E;
    pars->E_tot = f_Etot(pars);
    simp_v2(&theta_integrand_vec, Fnu, t, nu, Nt, 0.0, theta_wing, Ntheta,
                                    pars);
    double d_L = pars->d_L;
    double Fcoeff = cgs2mJy / (4*PI * d_L*d_L);
    
    int i;
    for(i=0; i<Nt; i++)
        Fnu[i] *= Fcoeff;
}

double flux_cone(double t_obs, double nu_obs, double E_iso, double theta_h,
                    double theta_cone_low, double theta_cone_hi, double atol,
                    struct fluxParams *pars)
{
    //printf("      t: %.3le th1: %.3f th2 %.3f\n", t_obs, theta_cone_low,
    //        theta_cone_hi);

    double theta_obs, theta_obs_cur;
    double F1, F2, Fboth;
    
    theta_obs = pars->theta_obs;

    if(E_iso > 0.0 && theta_h > 0.0)
        set_jet_params(pars, E_iso, theta_h);

    //Jet 
    theta_obs_cur = theta_obs;
    set_obs_params(pars, t_obs, nu_obs, theta_obs_cur, 
                    theta_cone_hi, theta_cone_low);
    F1 = flux(pars, atol);
    
    //Counter-jet
    theta_obs_cur = 180*deg2rad - theta_obs;
    set_obs_params(pars, t_obs, nu_obs, theta_obs_cur, 
                    theta_cone_hi, theta_cone_low);
    F2 = flux(pars, atol);

    Fboth = F1 + F2;

    if(F1 != F1 || F1 < 0.0)
        printf("bad F1:%.3lg t_obs=%.3le theta_lo=%.3lf theta_hi=%.3lf\n",
                F1, t_obs, theta_cone_low, theta_cone_hi);
    if(F2 != F2 || F2 < 0.0)
        printf("bad F2:%.3lg t_obs=%.3le theta_lo=%.3lf theta_hi=%.3lf\n",
                F2, t_obs, theta_cone_low, theta_cone_hi);

    return Fboth;
}

void calc_flux_density(int jet_type, int spec_type, double *t, double *nu,
                            double *Fnu, int N,
                            double theta_obs, double E_iso_core,
                            double theta_h_core, double theta_h_wing, 
                            double n_0, double p, double epsilon_E,
                            double epsilon_B, double ksi_N, double d_L,
                            int tRes, int latRes, double rtol, double *mask,
                            int nmask)
{
    double ta = t[0];
    double tb = t[0];
    int i;
    for(i=0; i<N; i++)
    {
        if(t[i] < ta)
            ta = t[i];
        else if(t[i] > tb)
            tb = t[i];
    }

    int res_cones = (int) (latRes*theta_h_wing / theta_h_core);

    struct fluxParams fp;
    setup_fluxParams(&fp, d_L, theta_obs, E_iso_core, theta_h_core,
                        theta_h_wing,
                        n_0, p, epsilon_E, epsilon_B, ksi_N, ta, tb, tRes,
                        spec_type, rtol, mask, nmask);

    if(jet_type == _tophat)
    {
        lc_tophat(t, nu, Fnu, N, E_iso_core, theta_h_core, &fp);
    }
    else if(jet_type == _cone)
    {
        lc_cone(t, nu, Fnu, N, E_iso_core, theta_h_core, theta_h_wing, &fp);
    }
    else if(jet_type == _powerlaw_core)
    {
        lc_powerlawCore(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, -2, NULL, NULL, res_cones, &fp);
    }
    else if(jet_type == _powerlaw)
    {
        lc_powerlaw(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &fp);
    }
    else if(jet_type == _Gaussian)
    {
        lc_Gaussian(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                    theta_h_wing, NULL, NULL, res_cones, &fp);
    }
    else if(jet_type == _Gaussian_core)
    {
        lc_GaussianCore(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                    theta_h_wing, NULL, NULL, res_cones, &fp);
    }
    else if(jet_type == _tophat + 10)
    {
        lc_vec(t, nu, Fnu, N, E_iso_core, theta_h_core, theta_h_core, 
                res_cones, &f_E_tophat, &f_Etot_tophat, &fp);
    }
    else if(jet_type == _Gaussian + 10)
    {
        lc_vec(t, nu, Fnu, N, E_iso_core, theta_h_core, theta_h_wing, 
                res_cones, &f_E_Gaussian, &f_Etot_Gaussian, &fp);
    }
    else if(jet_type == _powerlaw + 10)
    {
        lc_vec(t, nu, Fnu, N, E_iso_core, theta_h_core, theta_h_wing, 
                res_cones, &f_E_powerlaw, &f_Etot_powerlaw, &fp);
    }
    free_fluxParams(&fp);
}
///////////////////////////////////////////////////////////////////////////////

void setup_fluxParams(struct fluxParams *pars,
                        double d_L, double theta_obs, double E_iso_core,
                        double theta_core, double theta_wing, 
                        double n_0, double p, double epsilon_E,
                        double epsilon_B, double ksi_N, double ta, double tb,
                        double tRes, int spec_type, double flux_rtol,
                        double *mask, int nmask)
{
    pars->t_table = NULL;
    pars->R_table = NULL;
    pars->mu_table = NULL;
    pars->spec_type = spec_type;

    pars->d_L = d_L;
    pars->theta_obs = theta_obs;
    pars->E_iso_core = E_iso_core;
    pars->theta_core = theta_core;
    pars->theta_wing = theta_wing;
    pars->E_tot = -1;
    
    pars->n_0 = n_0;
    pars->p = p;
    pars->epsilon_E = epsilon_E;
    pars->epsilon_B = epsilon_B;
    pars->ksi_N = ksi_N;

    pars->ta = ta;
    pars->tb = tb;
    pars->tRes = tRes;
    pars->flux_rtol = flux_rtol;

    pars->mask = mask;
    pars->nmask = nmask;
}

///////////////////////////////////////////////////////////////////////////////

void set_jet_params(struct fluxParams *pars, double E_iso, double theta_h)
{
    double E_jet;
    if(pars->E_tot > 0.0)
        E_jet = pars->E_tot;
    else
        E_jet = (1.0 - cos(theta_h)) * E_iso;
    double n_0 = pars->n_0;
    double C_BM = sqrt(17.0 * E_iso / (8.0 * PI * m_p * n_0
                                        * pow( v_light, 5.0)));
    double C_ST = 2.0 / 5.0 * 1.15 * pow(E_jet / (m_p * n_0), 1.0 / 5.0 )
                            * invv_light;

    pars->E_iso = E_iso;
    pars->theta_h = theta_h;
    pars->C_BMsqrd = C_BM * C_BM;
    pars->C_STsqrd = C_ST * C_ST;
    pars->t_NR = pow(2.0, 1.0 / 3.0) * pow(C_BM, 2.0 / 3.0);

    make_R_table(pars);
}

///////////////////////////////////////////////////////////////////////////////

void set_obs_params(struct fluxParams *pars, double t_obs, double nu_obs,
                        double theta_obs_cur, double current_theta_cone_hi, 
                        double current_theta_cone_low)
{
    pars->t_obs = t_obs;
    pars->nu_obs = nu_obs;
    pars->theta_obs_cur = theta_obs_cur;
    pars->cto = cos(theta_obs_cur);
    pars->sto = sin(theta_obs_cur);
    pars->current_theta_cone_hi = current_theta_cone_hi;
    pars->current_theta_cone_low = current_theta_cone_low;
}

///////////////////////////////////////////////////////////////////////////////

void free_fluxParams(struct fluxParams *pars)
{
    if(pars->t_table != NULL)
    {
        free(pars->t_table);
        pars->t_table = NULL;
    }
    if(pars->R_table != NULL)
    {
        free(pars->R_table);
        pars->R_table = NULL;
    }
    if(pars->mu_table != NULL)
    {
        free(pars->mu_table);
        pars->mu_table = NULL;
    }
}

