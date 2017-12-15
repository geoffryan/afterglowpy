#include "offaxis_struct.h"

double dmin(const double a, const double b)
{
    if(a <= b)
        return a;
    else
        return b;
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

double get_t_e(double a_mu, double t_obs, double *mu_table, double *t_table,
                int table_entries)
{
  if(a_mu > mu_table[table_entries-1])
  {
    printf("mu >> 1? this should not have happened\n");
    printf("   t_obs=%.6lg mu=%.6lg mu_table[-1]=%.6lg, t_table[-1]=%.6lg\n",
                t_obs, a_mu, mu_table[table_entries-1], 
                t_table[table_entries-1]);
    abort();
  }

  if (mu_table[0] >= a_mu) // happens only if t_e very small
    return t_obs / (1.0 - a_mu); // so return small t_e limit
 
  // otherwise return linear interpolation between table entries
  unsigned int i = ((unsigned int) table_entries) >> 1;
  unsigned int a = 0;
  unsigned int b = table_entries-1;
  
  while (b-a > 1u)
  {
    i = (b+a) >> 1;
    if (mu_table[i] > a_mu)
        b = i;
    else
        a = i;
  }
  
  return ((a_mu - mu_table[a]) * t_table[b] + (mu_table[b] - a_mu) *
        t_table[a]) / (mu_table[b] - mu_table[a]);
}

///////////////////////////////////////////////////////////////////////////////

double Rintegrand(double a_t_e, void* params)
{
    double C_BMsqrd = ((double *) params)[0];
    double C_STsqrd = ((double *) params)[1];
    double lfbsqrd = get_lfacbetashocksqrd(a_t_e, C_BMsqrd, C_STsqrd);
    return v_light * sqrt(1.0 - 1.0 / (1.0 + lfbsqrd));
}

///////////////////////////////////////////////////////////////////////////////

double get_R(double a_t_e, double Rt0, double Rt1, double *R_table, 
                double *t_table, double *alpha_table, int table_entries)
{
// get R from table
  int i;

  // approximate and extrapolate outside of tabulated domain  
  if (a_t_e < t_table[2]) { return v_light * a_t_e; }
  if (a_t_e > t_table[table_entries - 3]) 
    return R_table[table_entries - 1] * pow(a_t_e / t_table[table_entries - 1], 
      2.0/5.0);

  // find nearest entry and return value assuming power law behaviour
  i = (int) ((table_entries - 1.0) * log(a_t_e / Rt0) / log(Rt1 / Rt0));
  return R_table[i] * pow(a_t_e / t_table[i], alpha_table[i]);
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
    double Rt0 = pars->Rt0;
    double Rt1 = pars->Rt1;
    double C_BMsqrd = pars->C_BMsqrd;
    double C_STsqrd = pars->C_STsqrd;
    double *t_table = pars->t_table;
    double *R_table = pars->R_table;
    double *alpha_table = pars->alpha_table;
    int table_entries = pars->table_entries;

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

    // set power law slopes at table times
    for (i=1; i < table_entries - 1; i++)
    {
        alpha_table[i] = 0.5 * (log(R_table[i] / R_table[i-1])
                                    / log(t_table[i] / t_table[i-1])
                                + log(R_table[i+1] / R_table[i])
                                    / log(t_table[i+1] / t_table[i]));
    }

    // free memory for integration routine
#ifdef USEGSL
    gsl_integration_workspace_free(w);
#endif
}

///////////////////////////////////////////////////////////////////////////////

double theta_integrand(double a_theta, void* params) // inner integral
{
    struct fluxParams *pars = (struct fluxParams *) params;
    
    //double cp = cos(pars->phi); 
    //double cto = cos(pars->theta_obs_cur);
    //double sto = sin(pars->theta_obs_cur);
    double ast = sin(a_theta);
    double act = cos(a_theta);

    double mu = ast * (pars->cp) * (pars->sto) + act * (pars->cto);

  double t_e = get_t_e(mu, pars->t_obs, pars->mu_table, pars->t_table, 
                        pars->table_entries);
  double R = get_R(t_e, pars->Rt0, pars->Rt1, pars->R_table, pars->t_table,
                        pars->alpha_table, pars->table_entries);
  //printf("%e, %e, %e # tobs, R, t_e\n", t_obs, t_e, R);
  double lfacbetashocksqrd = get_lfacbetashocksqrd(t_e, pars->C_BMsqrd, 
                                                    pars->C_STsqrd);
  double lfacbetasqrd = get_lfacbetasqrd(t_e, pars->C_BMsqrd, pars->C_STsqrd);
  
  // set remaining fluid quantities
  double lfacshocksqrd = 1.0 + lfacbetashocksqrd;
  double betashocksqrd = 1.0 - 1.0 / lfacshocksqrd;
  double betashock = sqrt(betashocksqrd);
  double lfacsqrd = 1.0 + lfacbetasqrd;
  double betasqrd = 1.0 - 1.0 / lfacsqrd;
  double lfac = sqrt( lfacsqrd );
  double beta = sqrt( betasqrd );
  double nprime = 4.0 * pars->n_0 * lfac; // comoving number density
  double e_th = (lfac - 1.0) * nprime * m_p * v_light * v_light;
  double B = sqrt(pars->epsilon_B * 8.0 * PI * e_th);
  double a = (1.0 - mu * beta); // beaming factor
  double ashock = (1.0 - mu * betashock); // shock velocity beaming factor
  double DR = R / (12.0 * lfacsqrd * ashock);
  if (DR < 0.0) DR *= -1.0; // DR is function of the absolute value of mu

  // set local emissivity 
  double p = pars->p;
  double nuprime = pars->nu_obs * lfac * a; // comoving observer frequency
  double lfac_m = (2.0 - p) / (1.0 - p) * (pars->epsilon_E * e_th / (
                    pars->ksi_N * nprime * m_e * v_light * v_light));
  double lfac_c = 6 * PI * m_e * lfac * v_light / (sigma_T * B * B * t_e);
  double nu_m = 3.0 * lfac_m * lfac_m * e_e * B / (4.0 * PI * m_e * v_light);
  double nu_c = 3.0 * lfac_c * lfac_c * e_e * B / (4.0 * PI * m_e * v_light);
  double em = pars->ksi_N * nprime * B;
  double freq = 0.0; // frequency dependent part of emissivity

  // set frequency dependence
  if (nu_c > nu_m)
  {
    if (nuprime < nu_m) 
      //freq = cbrt(nuprime / nu_m);
      freq = pow(nuprime / nu_m, 1.0 / 3.0 );
    else if (nuprime < nu_c)
      freq = pow(nuprime / nu_m, 0.5 * (1.0 - p));
    else
      freq = pow(nu_c / nu_m, 0.5 * (1.0 - p)) * pow(nuprime / nu_c, -0.5*p);
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

  return R * R * ast * DR * em * freq / (lfacsqrd * a * a);
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
  result = romb(&theta_integrand, theta_1-Dtheta, theta_1, 1000, 0, THETA_ACC,
                    params);
#endif
  
  //return result
  return result;
}

///////////////////////////////////////////////////////////////////////////////

double flux(struct fluxParams *pars) // determine flux for a given t_obs
{
  double result;
  double phi_0 = 0.0;
  double phi_1 = 2.0 * PI;
  
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
  
  //printf("about to integrate phi between %e and %e\n", phi_0, phi_1); fflush(stdout);
#ifdef USEGSL
  gsl_integration_qags (&F, phi_0, phi_1, 0, 1.0e-3, 1000, w, 
                            &result, &error); 
  // free memory
  gsl_integration_workspace_free(w);
#else
  result = romb(&phi_integrand, phi_0, phi_1, 1000, 0, PHI_ACC, pars);
#endif

  //return result
  double d_L = pars->d_L;
  double p = pars->p;

  return result / (4.0 * PI * d_L * d_L) * (p - 1.0) / 2.0 * sqrt(3.0) *
    e_e * e_e * e_e / (m_e * v_light * v_light) * cgs2mJy;
}

void lc_tophat(double *t, double *nu, double *F, int Nt,
                double E_iso, double theta_h, struct fluxParams *pars)
{
    int i;

    set_jet_params(pars, E_iso, theta_h);

    for(i=0; i<Nt; i++)
        F[i] = flux_cone(t[i], nu[i], -1, -1, 0.0, theta_h, pars);
}

void lc_powerlaw(double *t, double *nu, double *F, int Nt,
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
                                theta_cone_hi, pars);
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

    double Dtheta, theta_cone_hi, theta_cone_low, theta_h, theta_c, E_iso;

    if(theta_h_wing > 10*theta_h_core)
        theta_h_wing = 10*theta_h_core;

    Dtheta = theta_h_wing / res_cones;

    for(i=0; i<res_cones; i++)
    {
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
                                theta_cone_hi, pars);
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

    if(theta_h_wing > 10*theta_h_core)
        theta_h_wing = 10*theta_h_core;

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
                                theta_cone_hi, pars);
    }
}

double flux_cone(double t_obs, double nu_obs, double E_iso, double theta_h,
                    double theta_cone_low, double theta_cone_hi,
                    struct fluxParams *pars)
{
    double theta_obs, theta_obs_cur, theta_hi, theta_low;
    double F1, F2, Fboth;
    
    theta_obs = pars->theta_obs;
    theta_hi = theta_h;
    theta_low = 0.0;

    if(E_iso > 0.0 && theta_h > 0.0)
        set_jet_params(pars, E_iso, theta_h);

    //Jet 
    theta_obs_cur = theta_obs;
    set_obs_params(pars, t_obs, nu_obs, theta_obs_cur, 
                    theta_cone_hi, theta_cone_low);
    F1 = flux(pars);
    
    //Counter-jet
    theta_obs_cur = 180*deg2rad - theta_obs;
    set_obs_params(pars, t_obs, nu_obs, theta_obs_cur, 
                    theta_cone_hi, theta_cone_low);
    F2 = flux(pars);

    Fboth = F1 + F2;

    return Fboth;
}

void calc_flux_density(int jet_type, double *t, double *nu, double *Fnu, int N,
                            double theta_obs, double E_iso_core,
                            double theta_h_core, double theta_h_wing, 
                            double n_0, double p, double epsilon_E,
                            double epsilon_B, double ksi_N, double d_L)
{
    double Rt0 = 1.0e-2 * day2sec;
    double Rt1 = 1.0e7 * day2sec;
    int table_entries = 12000;
    int res_cones = 20;

    struct fluxParams fp;
    setup_fluxParams(&fp, d_L, theta_obs, n_0, p, epsilon_E, epsilon_B,
                        ksi_N, Rt0, Rt1, table_entries);

    if(jet_type == _tophat)
    {
        lc_tophat(t, nu, Fnu, N, E_iso_core, theta_h_core, &fp);
    }
    else if(jet_type == _powerlaw)
    {
        lc_powerlaw(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, -2, NULL, NULL, res_cones, &fp);
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
    free_fluxParams(&fp);
}
///////////////////////////////////////////////////////////////////////////////

void setup_fluxParams(struct fluxParams *pars, double d_L, double theta_obs,
                        double n_0, double p, double epsilon_E,
                        double epsilon_B, double ksi_N, double Rt0, double Rt1,
                        int table_entries)
{
    pars->t_table = (double *)malloc(sizeof(double) * table_entries);
    pars->R_table = (double *)malloc(sizeof(double) * table_entries);
    pars->mu_table = (double *)malloc(sizeof(double) * table_entries);
    pars->alpha_table = (double *)malloc(sizeof(double) * table_entries);
    pars->table_entries = table_entries;

    pars->d_L = d_L;
    pars->theta_obs = theta_obs;
    pars->n_0 = n_0;
    pars->p = p;
    pars->epsilon_E = epsilon_E;
    pars->epsilon_B = epsilon_B;
    pars->ksi_N = ksi_N;

    pars->Rt0 = Rt0;
    pars->Rt1 = Rt1;
}

///////////////////////////////////////////////////////////////////////////////

void set_jet_params(struct fluxParams *pars, double E_iso, double theta_h)
{
    double E_jet = theta_h * theta_h * E_iso / 2.0;
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
    free(pars->t_table);
    free(pars->R_table);
    free(pars->mu_table);
    free(pars->alpha_table);
}

