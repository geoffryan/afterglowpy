#include "offaxis_struct.h"
#include "shockEvolution.h"

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
        double b = pars->b;
        return pars->E_iso_core / pow(sqrt(1 + x*x/b), b);
    }
    return 0.0;
}

double f_E_GaussianCore(double theta, void *params)
{
    struct fluxParams *pars = (struct fluxParams *)params;
    if(theta <= pars->theta_wing)
    {
        double x = theta / pars->theta_core;
        return pars->E_iso_core * exp(-0.5*(x*x-1));
    }
    return 0.0;
}

double f_E_powerlawCore(double theta, void *params)
{
    struct fluxParams *pars = (struct fluxParams *)params;
    if(theta <= pars->theta_wing)
    {
        double x = theta / pars->theta_core;
        double b = pars->b;
        return pars->E_iso_core / pow(x, b);
    }
    return 0.0;
}

double f_E_twocomponent(double theta, void *params)
{
    struct fluxParams *pars = (struct fluxParams *)params;

    double th0 = 0.15;
    double th1 = 0.02;
    double b0 = 42;
    double b1 = 15;
    double e0 = 1.0;
    double e1 = 0.3;
    double f0 = 1.0;
    double f1 = 1.0;
    if(theta > th0)
        f0 = exp(-b0*(theta-th0));
    if(theta > th1)
        f1 = exp(-b1*(theta-th1));
    if(theta <= pars->theta_wing)
    {
        double x = theta / pars->theta_core;
        return pars->E_iso_core / (1 + x*x);
    }

    return pars->E_iso_core * (f0 + f1) / (e0 + e1);
}

double f_E_exponential(double theta, void *params)
{
    struct fluxParams *pars = (struct fluxParams *)params;
    if(theta <= pars->theta_wing)
    {
        double x = theta / pars->theta_core;
        double b = pars->b;
        return pars->E_iso_core * exp(-pow(x, b));
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
        return -1;
    }

    if(mu_table[0] >= mu) // happens only if t_e very small
    {
        printf("very small mu: mu=%.3lg, mu[0]=%.3lg\n", mu, mu_table[0]);
        //return t_obs / (1.0 - mu); // so return small t_e limit
        return -1;
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
    double *t_table_inner = pars->t_table_inner;
    double *R_table_inner = pars->R_table_inner;
    double *mu_table_inner = pars->mu_table_inner;
    int table_entries_inner = pars->table_entries_inner;

    int i;
    for (i = 0; i< table_entries; i++)
    {
        mu_table[i] = (t_table[i] - t_obs) / R_table[i] * v_light;
    }
    for (i = 0; i< table_entries_inner; i++)
    {
        mu_table_inner[i] = (t_table_inner[i] - t_obs)
                            / R_table_inner[i] * v_light;
    }
}

void make_R_tableInterp(struct fluxParams *pars)
{
    int tRes = pars->tRes;
    double Rt0 = pars->Rt0;
    double Rt1 = pars->Rt1;
    int table_entries = (int)(tRes * log10(Rt1/Rt0));
    
    pars->table_entries = table_entries;

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
    double Rpar[2] = {pars->C_BMsqrd, pars->C_STsqrd};
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

void make_R_table(struct fluxParams *pars)
{
    int tRes = pars->tRes;
    double Rt0 = pars->Rt0;
    double Rt1 = pars->Rt1;
    int table_entries = (int)(tRes * log10(Rt1/Rt0));
    double *temp;

    pars->table_entries_inner = pars->table_entries;
    pars->table_entries = table_entries;

    temp = pars->t_table_inner;
    pars->t_table_inner = pars->t_table;
    pars->t_table = (double *)realloc(temp, sizeof(double) * table_entries);
    temp = pars->R_table_inner;
    pars->R_table_inner = pars->R_table;
    pars->R_table = (double *)realloc(temp, sizeof(double) * table_entries);
    temp = pars->u_table_inner;
    pars->u_table_inner = pars->u_table;
    pars->u_table = (double *)realloc(temp, sizeof(double) * table_entries);
    temp = pars->th_table_inner;
    pars->th_table_inner = pars->th_table;
    pars->th_table = (double *)realloc(temp, sizeof(double) * table_entries);
    temp = pars->mu_table_inner;
    pars->mu_table_inner = pars->mu_table;
    pars->mu_table = (double *)realloc(temp, sizeof(double) * table_entries);

    double *t_table = pars->t_table;
    double *R_table = pars->R_table;
    double *u_table = pars->u_table;
    double *th_table = pars->th_table;

    double fac = pow(Rt1/Rt0, 1.0/(table_entries-1.0));
    t_table[0] = Rt0;

    int i;
    for(i=1; i<table_entries; i++)
        t_table[i] = t_table[i-1] * fac;

    double Rpar[2] = {pars->C_BMsqrd, pars->C_STsqrd};
    double R0 = romb(&Rintegrand, 0.0, Rt0, 1000, 0, R_ACC, Rpar);
    double u0 = sqrt(get_lfacbetasqrd(Rt0, pars->C_BMsqrd, pars->C_STsqrd));
    double th0 = pars->theta_h;
    double fom = 2*sin(0.5*th0)*sin(0.5*th0); //Fraction of solid angle in jet.

    double Mej_sph;
    if(pars->g_core > 1.0)
        Mej_sph = pars->E_iso / ((pars->g_init - 1.0) * v_light*v_light);
    else
        Mej_sph = 0.0;

    double thC = pars->theta_core;
    if(thC <= 0.0)
        thC = pars->theta_wing;
    //thC = 0.08; //th0;
    double thCg = pars->theta_core_global;
    if(thCg <= 0.0)
        thCg = thC;

    double args[12] = {pars->E_iso, Mej_sph, m_p*pars->n_0, 0.0, 0.0, 0.0, 
                        pars->L0, pars->q, pars->ts, thC, th0, thCg};
    int spread = pars->spread;
    //printf("t0=%.6le R0=%.6le u0=%.6le\n", Rt0, R0, u0);
    //shockInitDecel(Rt0, &R0, &u0, args);
    shockInitFind(Rt0, &R0, &u0, pars->tRes/10, args);
    //printf("t0=%.6le R0=%.6le u0=%.6le\n", Rt0, R0, u0);

    args[0] = pars->E_iso * fom;
    args[1] = Mej_sph * fom;
    shockEvolveSpreadRK4(t_table, R_table, u_table, th_table, table_entries,
                            R0, u0, th0, args, spread);

    if(R_table[0] != R_table[0])
    {
        printf("Rintegration Error: R[0]=%.3e  (fac=%.3e)\n", 
                                    R_table[0], fac);
        printf("    t=%.3e Rdot=%.3e t=%.3e Rdot=%.3e\n",
                    0.0,Rintegrand(0.0,Rpar), t_table[0], 
                    Rintegrand(t_table[0],Rpar));
    }
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
        printf("bad em:%.3le te=%.3le sinTheta=%.3lf mu=%.3lf\n",
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

    if(t_e < 0.0)
    {
        printf("BAD t_e: %.6lf Eiso=%.3le n0=%.3le thetah=%.3le\n",
                t_e, pars->E_iso, pars->n_0, pars->theta_h);
        printf("    theta_obs=%.3lf phi=%.3lf theta=%.3lf mu=%.3lf\n",
                pars->theta_obs, pars->phi, pars->theta, mu);
        printf("    L0=%.3le q=%.3lf ts=%.3le\n", pars->L0, pars->q, pars->ts);
        printf("    t[0]=%.3le t[-1]=%.3le R[0]=%.3le R[-1]=%.3le\n",
                pars->t_table[0], pars->t_table[pars->table_entries-1],
                pars->R_table[0], pars->R_table[pars->table_entries-1]);
        printf("    u[0]=%.3le u[-1]=%.3le th[0]=%.3le th[-1]=%.3le\n",
                pars->u_table[0], pars->u_table[pars->table_entries-1],
                pars->th_table[0], pars->th_table[pars->table_entries-1]);
        abort();
    }
    
    double R = interpolateLog(ia, ib, t_e, pars->t_table, pars->R_table, 
                            pars->table_entries);

    double us, u;
    if(pars->u_table != NULL)
    {
        u = interpolateLog(ia, ib, t_e, pars->t_table, pars->u_table,
                                        pars->table_entries);
        us = shockVel(u);
    }
    else
    {
        double us2 = get_lfacbetashocksqrd(t_e, pars->C_BMsqrd, pars->C_STsqrd);
        double u2 = get_lfacbetasqrd(t_e, pars->C_BMsqrd, pars->C_STsqrd);
        us = sqrt(us2);
        u = sqrt(u2);
    }
    
    double dFnu =  emissivity(pars->nu_obs, R, ast, mu, t_e, u, us,
                                pars->n_0, pars->p, pars->epsilon_E,
                                pars->epsilon_B, pars->ksi_N, pars->spec_type);

    if(dFnu != dFnu || dFnu < 0.0)
    {
        printf("bad dFnu:%.3le nu=%.3le R=%.3le th=%.3lf mu=%.3lf\n",
                dFnu, pars->nu_obs, R, a_theta, mu);
        printf("               t=%.3le u=%.3le us=%.3le n0=%.3le p=%.3lf\n",
                t_e, u, us, pars->n_0, pars->p);
        printf("               epse=%.3le epsB=%.3le ksiN=%.3le specType=%d\n",
                pars->epsilon_E, pars->epsilon_B, pars->ksi_N, pars->spec_type);
        printf("               Rt0=%.3le Rt1=%.3le E_iso=%.3le L0=%.3le ts=%.3le\n",
                pars->Rt0, pars->Rt1, pars->E_iso, pars->L0, pars->ts);
    }

    int i;
    double fac = 1.0;
    for(i=0; i<pars->nmask; i++)
    {
        double *m = &((pars->mask)[9*i]);
        if(m[0]<t_e && t_e<m[1] && m[2]<R && R<m[3] && m[4]<a_theta
                && a_theta<m[5] && m[6]<pars->phi && pars->phi<m[7])
            fac = m[8];
    }

    if(fac != fac || fac < 0.0)
        printf("bad mask fac: %.3le\n", fac);

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
    int spreadVersion = 1;
    if(pars->th_table != NULL && spreadVersion==1)
    {
        double th_0, th_1;
        th_1 = find_jet_edge(a_phi, pars->cto, pars->sto, theta_1,
                             pars->mu_table, pars->th_table,
                             pars->table_entries);

        if(0 || pars->table_entries_inner == 0)
        {
            double frac = theta_0 / theta_1;
            th_0 = frac * th_1;
        }
        else
        {
            th_0 = find_jet_edge(a_phi, pars->cto, pars->sto, theta_0,
                                 pars->mu_table_inner, pars->th_table_inner,
                                 pars->table_entries_inner);
        }
        /*
        double frac = theta_0 / theta_1;
        double th_0 = frac * th_1;
        */
        theta_0 = th_0;
        theta_1 = th_1;
        if(theta_0 > 0.5*M_PI)
            theta_0 = 0.5*M_PI;
        if(theta_1 > 0.5*M_PI)
            theta_1 = 0.5*M_PI;
    }
    if(pars->th_table != NULL && spreadVersion==2)
    {
        // approx mu
        double ct = cos(0.5*(theta_0+theta_1));
        double st = sin(0.5*(theta_0+theta_1));
        double mu = pars->cp * st * (pars->sto) + ct * (pars->cto);

        int ia = searchSorted(mu, pars->mu_table, pars->table_entries);
        int ib = ia+1;
        double th = interpolateLin(ia, ib, mu, pars->mu_table, pars->th_table, 
                                    pars->table_entries);

        theta_0 *= th/pars->theta_h;
        theta_1 *= th/pars->theta_h;
        if(theta_0 > 0.5*M_PI)
            theta_0 = 0.5*M_PI;
        if(theta_1 > 0.5*M_PI)
            theta_1 = 0.5*M_PI;
    }
    else if (pars->t_obs > pars->t_NR && spreadVersion==3)
    {
        theta_1 = dmin(0.5 * PI, 
                        pars->theta_h + 0.1 * log( pars->t_obs / pars->t_NR));
        if(theta_0 != 0.0)
            theta_0 = theta_1-Dtheta;
    }

    if(theta_0 >= theta_1)
        return 0.0;
 
    //printf("# theta integration domain: %e - %e\n", theta_1 - Dtheta, theta_1); fflush(stdout);
 
    // For a given phi, integrate over theta
#ifdef USEGSL
    gsl_integration_qags(&F, theta_0, theta_1, 0, 1.0e-4, 1000, w, 
                            &result, &error);
  // free integration routine memory
    gsl_integration_workspace_free(w);
#else
    result = romb(&theta_integrand, theta_0, theta_1, 1000, 
                        pars->theta_atol, THETA_ACC, params);
#endif
    if(result != result || result < 0.0)
        printf("bad result:%.3le t_obs=%.3le theta_lo=%.3lf theta_hi=%.3lf theta_log=%.3lf phi=%.3lf\n",
              result, pars->t_obs, theta_0, theta_1, pars->theta_h + 0.1 * log( pars->t_obs / pars->t_NR), pars->phi);
  
    //return result
    return result;
}

double find_jet_edge(double phi, double cto, double sto, double theta0,
                     double *a_mu, double *a_thj, int N)
{
    double cp = cos(phi);
    double mu = cos(theta0)*cto + sin(theta0)*sto*cp;

    int ia = searchSorted(mu, a_mu, N);
    
    if(a_thj[ia] <= theta0 && theta0 <= a_thj[ia+1])
        return theta0;

    double tha, thb;
    if(theta0 < a_thj[ia])
    {
        //The jet is spreading
        tha = theta0;
        thb = 0.5*M_PI;
    }
    else
    {
        //Guessed too far out!
        tha = 0.0;
        thb = theta0;
    }

    int i = 0;
    while(thb-tha > 1.0e-5 && i < 100)
    {
        double th = 0.5*(tha+thb);
        mu = cos(th)*cto + sin(th)*sto*cp;
        ia = searchSorted(mu, a_mu, N);
        if(th < a_thj[ia])
            tha = th;
        else
            thb = th;
        i++;
    }

    //printf("iter: %d, th0=%.6lf, th=%.6lf\n", i, theta0, tha);

    return tha;
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
  //result = 2 * Fcoeff * PI * phi_integrand(0.0, pars);
#endif

  //return result

  return result;
}

void lc_cone(double *t, double *nu, double *F, int Nt, double E_iso,
                double theta_core, double theta_wing, struct fluxParams *pars)
{
    int i;

    set_jet_params(pars, E_iso, theta_wing);

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

void lc_struct(double *t, double *nu, double *F, int Nt,
                        double E_iso_core, 
                        double theta_h_core, double theta_h_wing,
                        double *theta_c_arr, double *E_iso_arr,
                        int res_cones, double (*f_E)(double,void *),
                        struct fluxParams *pars)
{
    //Flux from a structured jet.
    
    int i,j;
    //No Core
    for(j=0; j<Nt; j++)
        F[j] = 0.0;

    double Dtheta, theta_cone_hi, theta_cone_low, theta_h, theta_c, E_iso;

    Dtheta = theta_h_wing / res_cones;

    for(i=0; i<res_cones; i++)
    {
        theta_c = (i+0.5) * Dtheta;
        E_iso = f_E(theta_c, pars);

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

void lc_structCore(double *t, double *nu, double *F, int Nt,
                        double E_iso_core, 
                        double theta_h_core, double theta_h_wing,
                        double *theta_c_arr, double *E_iso_arr,
                        int res_cones, double (*f_E)(double,void *),
                        struct fluxParams *pars)
{
    //Flux from a structured jet with core.
    
    lc_tophat(t, nu, F, Nt, E_iso_core, theta_h_core, pars);

    double Dtheta, theta_cone_hi, theta_cone_low, theta_h, theta_c, E_iso;

    Dtheta = (theta_h_wing - theta_h_core) / res_cones;

    int i,j;
    for(i=0; i<res_cones; i++)
    {
        theta_c = theta_h_core + (i+0.5) * Dtheta;
        E_iso = f_E(theta_c, pars);

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
    /*
    theta_obs_cur = 180*deg2rad - theta_obs;
    set_obs_params(pars, t_obs, nu_obs, theta_obs_cur, 
                    theta_cone_hi, theta_cone_low);
    F2 = flux(pars, atol);
    */
    F2 = 0.0;
    Fboth = F1 + F2;


    if(F1 != F1 || F1 < 0.0)
        printf("bad F1:%.3lg t_obs=%.3le theta_lo=%.3lf theta_hi=%.3lf\n",
                F1, t_obs, theta_cone_low, theta_cone_hi);
    if(F2 != F2 || F2 < 0.0)
        printf("bad F2:%.3lg t_obs=%.3le theta_lo=%.3lf theta_hi=%.3lf\n",
                F2, t_obs, theta_cone_low, theta_cone_hi);

    return Fboth;
}

double intensity(double theta, double phi, double tobs, double nuobs,
                double theta_obs, double theta_cone_hi, double theta_cone_low,
                struct fluxParams *pars)
{
    double I = 0;

    int remakeMu = 0;
    if(tobs != pars->t_obs)
        remakeMu = 1;
    set_obs_params(pars, tobs, nuobs, theta_obs, theta_cone_hi, 
                    theta_cone_low);

    if(remakeMu)
        make_mu_table(pars);

    double mu = cos(theta)*cos(theta_obs)
                    + sin(theta)*sin(theta_obs)*cos(phi);

    int ia = searchSorted(mu, pars->mu_table, pars->table_entries);
    int ib = ia + 1;
    double t_e = interpolateLin(ia, ib, mu, pars->mu_table,
                                pars->t_table, pars->table_entries);
    t_e = check_t_e(t_e, mu, pars->t_obs, pars->mu_table, 
                                pars->table_entries);
    if(t_e < 0.0)
        printf("WTFWTF\n");

    double R = interpolateLog(ia, ib, t_e, pars->t_table,
                                pars->R_table, pars->table_entries);
    double u = interpolateLog(ia, ib, t_e, pars->t_table,
                                pars->u_table, pars->table_entries);
    double us = shockVel(u);

    I = emissivity(pars->nu_obs, R, 1.0, mu, t_e, u, us, pars->n_0,
                        pars->p, pars->epsilon_E, pars->epsilon_B, 
                        pars->ksi_N, pars->spec_type);

    return I;
}

void shockVals(double theta, double phi, double tobs,
                 double *t, double *R, double *u, double *thj,
                 double theta_obs, double theta_cone_hi, double theta_cone_low,
                 struct fluxParams *pars)
{
    int remakeMu = 0;
    if(tobs != pars->t_obs)
        remakeMu = 1;
    set_obs_params(pars, tobs, 1.0, theta_obs, theta_cone_hi, 
                    theta_cone_low);

    if(remakeMu)
        make_mu_table(pars);

    double mu = cos(theta)*cos(theta_obs)
                    + sin(theta)*sin(theta_obs)*cos(phi);

    int ia = searchSorted(mu, pars->mu_table, pars->table_entries);
    int ib = ia + 1;
    double t_e = interpolateLin(ia, ib, mu, pars->mu_table,
                                pars->t_table, pars->table_entries);
    t_e = check_t_e(t_e, mu, pars->t_obs, pars->mu_table, 
                                pars->table_entries);
    if(t_e < 0.0)
        printf("WTFWTF\n");


    *t = t_e;
    *R = interpolateLog(ia, ib, t_e, pars->t_table,
                        pars->R_table, pars->table_entries);
    *u = interpolateLog(ia, ib, t_e, pars->t_table,
                        pars->u_table, pars->table_entries);
    *thj = interpolateLin(ia, ib, t_e, pars->t_table,
                          pars->th_table, pars->table_entries);
}

void intensity_cone(double *theta, double *phi, double *t, double *nu, 
                        double *I, int N, double E_iso_core, 
                        double theta_h_core, double theta_h_wing, 
                        struct fluxParams *pars)
{
    //Intensity of a cone segment.
    
    int j;
    for(j=0; j<N; j++)
        I[j] = 0.0;

    double theta_obs = pars->theta_obs;
    double dL = pars->d_L;

    double Fcoeff = cgs2mJy / (4*M_PI*dL*dL);

    double theta_cone_hi = theta_h_wing;
    double theta_cone_low = theta_h_core;

    set_jet_params(pars, E_iso_core, theta_h_wing);
    set_obs_params(pars, t[0], nu[0], theta_obs, theta_cone_hi, theta_cone_low);
    make_mu_table(pars);

    for(j=0; j<N; j++)
    {
        double th = theta[j];
        double ph = phi[j];

        if(I[j] > 0.0 || th < theta_cone_low)
            continue;

        double th_a, th_b;
        th_b = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                             pars->mu_table, pars->th_table,
                             pars->table_entries);
        if(pars->table_entries_inner == 0)
            th_a = (theta_cone_low / theta_cone_hi) * th_b;
        else
            th_a = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                 pars->mu_table_inner, 
                                 pars->th_table_inner,
                                 pars->table_entries_inner);

        if(th < th_a || th > th_b)
            continue;

        I[j] = Fcoeff * intensity(th, ph, t[j], nu[j], theta_obs,        
                                    theta_cone_hi, theta_cone_low, pars);
    }
}

void intensity_struct(double *theta, double *phi, double *t, double *nu, 
                        double *I, int N,
                        double E_iso_core, 
                        double theta_h_core, double theta_h_wing,
                        int res_cones, double (*f_E)(double,void *),
                        struct fluxParams *pars)
{
    //Intensity of a structured jet.
    
    int i,j;
    //No Core
    for(j=0; j<N; j++)
        I[j] = 0.0;

    double Dtheta, theta_cone_hi, theta_cone_low, theta_h, theta_c, E_iso;
    double theta_obs, dL;

    Dtheta = theta_h_wing / res_cones;
    theta_obs = pars->theta_obs;
    dL = pars->d_L;

    double Fcoeff = cgs2mJy / (4*M_PI*dL*dL);

    for(i=0; i<res_cones; i++)
    {
        theta_c = (i+0.5) * Dtheta;
        E_iso = f_E(theta_c, pars);

        theta_cone_hi = (i+1) * Dtheta;
        theta_cone_low = i * Dtheta;
        theta_h = theta_cone_hi;

        set_jet_params(pars, E_iso, theta_h);
        set_obs_params(pars, t[0], nu[0], theta_obs, theta_cone_hi, 
                        theta_cone_low);
        make_mu_table(pars);

        for(j=0; j<N; j++)
        {
            double th = theta[j];
            double ph = phi[j];

            if(I[j] > 0.0 || th < theta_cone_low)
                continue;

            double th_a, th_b;
            th_b = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                 pars->mu_table, pars->th_table,
                                 pars->table_entries);
            if(pars->table_entries_inner == 0)
                th_a = (theta_cone_low / theta_cone_hi) * th_b;
            else
                th_a = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                     pars->mu_table_inner, 
                                     pars->th_table_inner,
                                     pars->table_entries_inner);

            if(th < th_a || th > th_b)
                continue;

            I[j] += Fcoeff * intensity(th, ph, t[j], nu[j], theta_obs,        
                                    theta_cone_hi, theta_cone_low, pars);
        }
    }
}

void intensity_structCore(double *theta, double *phi, double *t, double *nu, 
                            double *I, int N,
                            double E_iso_core, 
                            double theta_h_core, double theta_h_wing,
                            int res_cones, double (*f_E)(double,void *),
                            struct fluxParams *pars)
{
    //Intensity of a structured jet.
    
    int i,j;
    //Core
    intensity_cone(theta, phi, t, nu, I, N, E_iso_core, 0.0, theta_h_core,
                    pars);

    double Dtheta, theta_cone_hi, theta_cone_low, theta_h, theta_c, E_iso;
    double theta_obs, dL;

    Dtheta = theta_h_wing / res_cones;
    theta_obs = pars->theta_obs;
    dL = pars->d_L;

    double Fcoeff = cgs2mJy / (4*M_PI*dL*dL);

    for(i=0; i<res_cones; i++)
    {
        theta_c = (i+0.5) * Dtheta;
        E_iso = f_E(theta_c, pars);

        theta_cone_hi = (i+1) * Dtheta;
        theta_cone_low = i * Dtheta;
        theta_h = theta_cone_hi;

        set_jet_params(pars, E_iso, theta_h);
        set_obs_params(pars, t[0], nu[0], theta_obs, theta_cone_hi, 
                        theta_cone_low);
        make_mu_table(pars);
        
        for(j=0; j<N; j++)
        {
            double th = theta[j];
            double ph = phi[j];
            
            if(I[j] > 0.0 || th < theta_cone_low)
                continue;

            double th_a, th_b;
            th_b = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                 pars->mu_table, pars->th_table,
                                 pars->table_entries);
            if(pars->table_entries_inner == 0)
                th_a = (theta_cone_low / theta_cone_hi) * th_b;
            else
                th_a = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                     pars->mu_table_inner, 
                                     pars->th_table_inner,
                                     pars->table_entries_inner);

            if(th < th_a || th > th_b)
                continue;

            I[j] += Fcoeff * intensity(th, ph, t[j], nu[j], theta_obs,        
                                    theta_cone_hi, theta_cone_low, pars);
        }
    }
}

void shockVals_cone(double *theta, double *phi, double *tobs,
                    double *t, double *R, double *u, double *thj, 
                    int N, double E_iso_core, 
                    double theta_h_core, double theta_h_wing, 
                    struct fluxParams *pars)
{
    //Intensity of a cone segment.
    
    int j;
    for(j=0; j<N; j++)
    {
        t[j] = 0.0;
        R[j] = 0.0;
        u[j] = 0.0;
        thj[j] = 0.0;
    }

    double theta_obs = pars->theta_obs;

    double theta_cone_hi = theta_h_wing;
    double theta_cone_low = theta_h_core;

    set_jet_params(pars, E_iso_core, theta_h_wing);
    set_obs_params(pars, tobs[0], 1.0, theta_obs,
                   theta_cone_hi, theta_cone_low);
    make_mu_table(pars);

    for(j=0; j<N; j++)
    {
        double th = theta[j];
        double ph = phi[j];

        if(t[j] > 0.0 || th < theta_cone_low)
            continue;

        double th_a, th_b;
        th_b = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                             pars->mu_table, pars->th_table,
                             pars->table_entries);
        if(pars->table_entries_inner == 0)
            th_a = (theta_cone_low / theta_cone_hi) * th_b;
        else
            th_a = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                 pars->mu_table_inner, 
                                 pars->th_table_inner,
                                 pars->table_entries_inner);

        if(th < th_a || th > th_b)
            continue;

        shockVals(th, ph, tobs[j], t+j, R+j, u+j, thj+j,
                    theta_obs, theta_cone_hi, theta_cone_low, pars);
    }
}

void shockVals_struct(double *theta, double *phi, double *tobs,
                        double *t, double *R, double *u, double *thj, int N,
                        double E_iso_core, 
                        double theta_h_core, double theta_h_wing,
                        int res_cones, double (*f_E)(double,void *),
                        struct fluxParams *pars)
{
    //Intensity of a structured jet.
    
    int i,j;
    //No Core
    for(j=0; j<N; j++)
    {
        t[j] = 0.0;
        R[j] = 0.0;
        u[j] = 0.0;
        thj[j] = 0.0;
    }

    double Dtheta, theta_cone_hi, theta_cone_low, theta_h, theta_c, E_iso;
    double theta_obs;

    Dtheta = theta_h_wing / res_cones;
    theta_obs = pars->theta_obs;

    for(i=0; i<res_cones; i++)
    {
        theta_c = (i+0.5) * Dtheta;
        E_iso = f_E(theta_c, pars);

        theta_cone_hi = (i+1) * Dtheta;
        theta_cone_low = i * Dtheta;
        theta_h = theta_cone_hi;

        set_jet_params(pars, E_iso, theta_h);
        set_obs_params(pars, tobs[0], 1.0, theta_obs, theta_cone_hi, 
                        theta_cone_low);
        make_mu_table(pars);

        for(j=0; j<N; j++)
        {
            double th = theta[j];
            double ph = phi[j];

            if(t[j] > 0.0 || th < theta_cone_low)
                continue;

            double th_a, th_b;
            th_b = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                 pars->mu_table, pars->th_table,
                                 pars->table_entries);
            if(pars->table_entries_inner == 0)
                th_a = (theta_cone_low / theta_cone_hi) * th_b;
            else
                th_a = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                     pars->mu_table_inner, 
                                     pars->th_table_inner,
                                     pars->table_entries_inner);

            if(th < th_a || th > th_b)
                continue;

            shockVals(th, ph, tobs[j], t+j, R+j, u+j, thj+j,
                      theta_obs, theta_cone_hi, theta_cone_low, pars);
        }
    }
}

void shockVals_structCore(double *theta, double *phi, double *tobs, 
                          double *t, double *R, double *u, double *thj, int N,
                          double E_iso_core, 
                          double theta_h_core, double theta_h_wing,
                          int res_cones, double (*f_E)(double,void *),
                          struct fluxParams *pars)
{
    //Intensity of a structured jet.
    
    int i,j;
    //Core
    shockVals_cone(theta, phi, tobs, t, R, u, thj, N, E_iso_core, 0.0,
                    theta_h_core, pars);

    double Dtheta, theta_cone_hi, theta_cone_low, theta_h, theta_c, E_iso;
    double theta_obs;

    Dtheta = theta_h_wing / res_cones;
    theta_obs = pars->theta_obs;

    for(i=0; i<res_cones; i++)
    {
        theta_c = (i+0.5) * Dtheta;
        E_iso = f_E(theta_c, pars);

        theta_cone_hi = (i+1) * Dtheta;
        theta_cone_low = i * Dtheta;
        theta_h = theta_cone_hi;

        set_jet_params(pars, E_iso, theta_h);
        set_obs_params(pars, tobs[0], 1.0, theta_obs, theta_cone_hi, 
                        theta_cone_low);
        make_mu_table(pars);
        
        for(j=0; j<N; j++)
        {
            double th = theta[j];
            double ph = phi[j];
            
            if(t[j] > 0.0 || th < theta_cone_low)
                continue;

            double th_a, th_b;
            th_b = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                 pars->mu_table, pars->th_table,
                                 pars->table_entries);
            if(pars->table_entries_inner == 0)
                th_a = (theta_cone_low / theta_cone_hi) * th_b;
            else
                th_a = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                     pars->mu_table_inner, 
                                     pars->th_table_inner,
                                     pars->table_entries_inner);

            if(th < th_a || th > th_b)
                continue;

            shockVals(th, ph, tobs[j], t+j, R+j, u+j, thj+j,
                      theta_obs, theta_cone_hi, theta_cone_low, pars);
        }
    }
}

void calc_flux_density(int jet_type, int spec_type, double *t, double *nu,
                            double *Fnu, int N,
                            double theta_obs, double E_iso_core,
                            double theta_h_core, double theta_h_wing, 
                            double b, double L0, double q, double ts, 
                            double n_0, double p, double epsilon_E,
                            double epsilon_B, double ksi_N, double d_L,
                            double g0, double E_core_global,
                            double theta_h_core_global,
                            int tRes, int latRes, double rtol, double *mask,
                            int nmask, int spread, int gamma_type)
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
                        theta_h_wing, b, L0, q, ts,
                        n_0, p, epsilon_E, epsilon_B, ksi_N, g0, 
                        E_core_global, theta_h_core_global, ta, tb, tRes,
                        spec_type, rtol, mask, nmask, spread, gamma_type);

    if(jet_type == _tophat)
    {
        lc_tophat(t, nu, Fnu, N, E_iso_core, theta_h_core, &fp);
    }
    else if(jet_type == _cone)
    {
        lc_cone(t, nu, Fnu, N, E_iso_core, theta_h_core, theta_h_wing, &fp);
    }
    else if(jet_type == _Gaussian)
    {
        lc_struct(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_Gaussian, &fp);
    }
    else if(jet_type == _powerlaw)
    {
        lc_struct(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_powerlaw, &fp);
    }
    else if(jet_type == _twocomponent)
    {
        lc_struct(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_twocomponent, &fp);
    }
    else if(jet_type == _Gaussian_core)
    {
        lc_structCore(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_GaussianCore, &fp);
    }
    else if(jet_type == _powerlaw_core)
    {
        lc_structCore(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_powerlawCore, &fp);
    }
    else if(jet_type == _exponential)
    {
        lc_structCore(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_exponential, &fp);
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

void calc_intensity(int jet_type, int spec_type, double *theta, double *phi,
                            double *t, double *nu, double *Inu, int N,
                            double theta_obs, double E_iso_core,
                            double theta_h_core, double theta_h_wing, 
                            double b, double L0, double q, double ts, 
                            double n_0, double p, double epsilon_E,
                            double epsilon_B, double ksi_N, double d_L,
                            double g0, double E_core_global,
                            double theta_h_core_global,
                            int tRes, int latRes, double rtol, double *mask,
                            int nmask, int spread, int gamma_type)
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
                        theta_h_wing, b, L0, q, ts,
                        n_0, p, epsilon_E, epsilon_B, ksi_N, g0, 
                        E_core_global, theta_h_core_global, ta, tb, tRes,
                        spec_type, rtol, mask, nmask, spread, gamma_type);

    if(jet_type == _tophat)
    {
        intensity_cone(theta, phi, t, nu, Inu, N, E_iso_core, 0.0, theta_h_core,
                            &fp);
    }
    else if(jet_type == _cone)
    {
        intensity_cone(theta, phi, t, nu, Inu, N, E_iso_core, theta_h_core,
                        theta_h_wing, &fp);
    }
    else if(jet_type == _Gaussian)
    {
        intensity_struct(theta, phi, t, nu, Inu, N, E_iso_core, theta_h_core, 
                theta_h_wing, res_cones, &f_E_Gaussian, &fp);
    }
    else if(jet_type == _powerlaw)
    {
        intensity_struct(theta, phi, t, nu, Inu, N, E_iso_core, theta_h_core, 
                theta_h_wing, res_cones, &f_E_powerlaw, &fp);
    }
    else if(jet_type == _Gaussian_core)
    {
        intensity_structCore(theta, phi, t, nu, Inu, N, E_iso_core, 
                            theta_h_core, theta_h_wing, res_cones,
                            &f_E_GaussianCore, &fp);
    }
    else if(jet_type == _powerlaw_core)
    {
        intensity_structCore(theta, phi, t, nu, Inu, N, E_iso_core, 
                                theta_h_core, theta_h_wing, res_cones, 
                                &f_E_powerlawCore, &fp);
    }
    free_fluxParams(&fp);

}

void calc_shockVals(int jet_type, double *theta, double *phi, double *tobs,
                    double *t, double *R, double *u, double *thj, int N,
                    double theta_obs, double E_iso_core,
                    double theta_h_core, double theta_h_wing, 
                    double b, double L0, double q, double ts, 
                    double n_0, double p, double epsilon_E,
                    double epsilon_B, double ksi_N, double d_L,
                    double g0, double E_core_global,
                    double theta_h_core_global,
                    int tRes, int latRes, double rtol, double *mask,
                    int nmask, int spread, int gamma_type)
{
    double ta = tobs[0];
    double tb = tobs[0];
    int i;
    for(i=0; i<N; i++)
    {
        if(tobs[i] < ta)
            ta = tobs[i];
        else if(tobs[i] > tb)
            tb = tobs[i];
    }

    int res_cones = (int) (latRes*theta_h_wing / theta_h_core);

    struct fluxParams fp;
    setup_fluxParams(&fp, d_L, theta_obs, E_iso_core, theta_h_core,
                        theta_h_wing, b, L0, q, ts,
                        n_0, p, epsilon_E, epsilon_B, ksi_N, g0, 
                        E_core_global, theta_h_core_global, ta, tb, tRes,
                        0, rtol, mask, nmask, spread, gamma_type);

    if(jet_type == _tophat)
    {
        shockVals_cone(theta, phi, tobs, t, R, u, thj, N, E_iso_core, 0.0,
                        theta_h_core, &fp);
    }
    else if(jet_type == _cone)
    {
        shockVals_cone(theta, phi, tobs, t, R, u, thj, N, E_iso_core,
                        theta_h_core, theta_h_wing, &fp);
    }
    else if(jet_type == _Gaussian)
    {
        shockVals_struct(theta, phi, tobs, t, R, u, thj, N, E_iso_core,
                         theta_h_core, theta_h_wing, res_cones, &f_E_Gaussian,
                         &fp);
    }
    else if(jet_type == _powerlaw)
    {
        shockVals_struct(theta, phi, tobs, t, R, u, thj, N, E_iso_core,
                         theta_h_core, theta_h_wing, res_cones, &f_E_powerlaw,
                         &fp);
    }
    else if(jet_type == _Gaussian_core)
    {
        shockVals_structCore(theta, phi, tobs, t, R, u, thj, N, E_iso_core,
                             theta_h_core, theta_h_wing, res_cones,
                             &f_E_GaussianCore, &fp);
    }
    else if(jet_type == _powerlaw_core)
    {
        shockVals_structCore(theta, phi, tobs, t, R, u, thj, N, E_iso_core,
                             theta_h_core, theta_h_wing, res_cones,
                             &f_E_powerlawCore, &fp);
    }
    free_fluxParams(&fp);

}
///////////////////////////////////////////////////////////////////////////////

void setup_fluxParams(struct fluxParams *pars,
                        double d_L, double theta_obs, 
                        double E_iso_core, 
                        double theta_core, double theta_wing, double b,
                        double L0, double q, double ts, 
                        double n_0, double p, double epsilon_E,
                        double epsilon_B, double ksi_N, double g0,
                        double E_core_global, double theta_core_global, 
                        double ta, double tb,
                        double tRes, int spec_type, double flux_rtol,
                        double *mask, int nmask, int spread, int gamma_type)
{
    pars->t_table = NULL;
    pars->R_table = NULL;
    pars->u_table = NULL;
    pars->th_table = NULL;
    pars->mu_table = NULL;
    pars->table_entries = 0;
    pars->t_table_inner = NULL;
    pars->R_table_inner = NULL;
    pars->u_table_inner = NULL;
    pars->th_table_inner = NULL;
    pars->mu_table_inner = NULL;
    pars->table_entries_inner = 0;

    pars->spec_type = spec_type;
    pars->gamma_type = gamma_type;

    pars->d_L = d_L;
    pars->theta_obs = theta_obs;
    pars->E_iso_core = E_iso_core;
    pars->theta_core = theta_core;
    pars->theta_wing = theta_wing;
    pars->b = b;
    pars->E_tot = -1;
    pars->g_core = g0;
    pars->E_core_global = E_core_global;
    pars->theta_core_global = theta_core_global;
   
    pars->L0 = L0;
    pars->q = q;
    pars->ts = ts;

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
    pars->spread = spread;
}

///////////////////////////////////////////////////////////////////////////////

void set_jet_params(struct fluxParams *pars, double E_iso, double theta_h)
{
    // min/max observer times
    double ta = pars->ta;
    double tb = pars->tb;

    double Rt0, Rt1;

    //at fixed t_obs, earliest emission is *always* from mu=-1
    // so t_obs ~ t_e
    Rt0 = 0.1*ta;
    double E_jet;
    if(pars->E_tot > 0.0)
        E_jet = pars->E_tot;
    else
        E_jet = (1.0 - cos(theta_h)) * E_iso;

    double Einj = 0.0;
    double ti = 0.0;
    if(pars->L0 > 0.0 && pars->ts > 0.0)
    {
        // Energy injection uses 1e3s as reference time. 
        Einj = E_inj(pars->ts, pars->L0, pars->q, pars->ts);
        ti = pars->ts;
    }
    E_jet += Einj;

    double c5 = v_light*v_light*v_light*v_light*v_light;
    double n_0 = pars->n_0;
    double C_BM = sqrt(17.0 * E_iso / (8.0 * PI * m_p * n_0 * c5));
    double C_ST = 2.0 / 5.0 * 1.15 * pow(E_jet / (m_p * n_0), 1.0 / 5.0 )
                            * invv_light;
    double t_NR = pow(2.0, 1.0 / 3.0) * pow(C_BM, 2.0 / 3.0);

    //if(Einj > 0.0)
    //    t_NR *= pow((E_iso+Einj)/E_iso, 1.0/3.0);
 

    pars->E_iso = E_iso;
    pars->theta_h = theta_h;
    if(pars->gamma_type == 1 && pars->E_core_global > 0.0)
        pars->g_init = 1.0 + (pars->g_core - 1) * E_iso / pars->E_core_global;
    else
        pars->g_init = 1.0 + (pars->g_core - 1) * E_iso / pars->E_iso_core;
    pars->C_BMsqrd = C_BM * C_BM;
    pars->C_STsqrd = C_ST * C_ST;
    pars->t_NR = t_NR;

    //This *should* be an over-estimate of the non-relativistic time.
    double t_NR2 = pow((E_iso+Einj) / (m_p*n_0 * c5), 1.0/3.0); 

    /*
    if(pars->L0 > 0.0 && pars->ts > 0.0)
    {
        if(Rt0 * pars->L0 > 0.1*E_iso)
            Rt0 = 0.1 * E_iso / pars->L0;

        C_BM *= sqrt((E_iso+Einj)/E_iso);
    }
    */

    //at fixed t_obs, latest emission is *always* from mu=+1
    // so t_obs ~ t-R/c
    /*
    if(tb > 0.1*t_NR)  // at late times R ~ c t_NR << c t so t_obs ~ t_e
        Rt1 = 100*(tb+t_NR);
    else // at early times t_obs ~ t*(gamma_sh^-2)/8 ~ CBM^-2 * t^4 / 8
        Rt1 = 100*pow(8*tb*C_BM*C_BM, 0.25);
    */

    Rt1 = 100*(tb + t_NR2 + ti);

    pars->Rt0 = Rt0;
    pars->Rt1 = Rt1;
    
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
    if(pars->u_table != NULL)
    {
        free(pars->u_table);
        pars->u_table = NULL;
    }
    if(pars->th_table != NULL)
    {
        free(pars->th_table);
        pars->th_table = NULL;
    }
    if(pars->mu_table != NULL)
    {
        free(pars->mu_table);
        pars->mu_table = NULL;
    }

    if(pars->t_table_inner != NULL)
    {
        free(pars->t_table_inner);
        pars->t_table_inner = NULL;
    }
    if(pars->R_table_inner != NULL)
    {
        free(pars->R_table_inner);
        pars->R_table_inner = NULL;
    }
    if(pars->u_table_inner != NULL)
    {
        free(pars->u_table_inner);
        pars->u_table_inner = NULL;
    }
    if(pars->th_table_inner != NULL)
    {
        free(pars->th_table_inner);
        pars->th_table_inner = NULL;
    }
    if(pars->mu_table_inner != NULL)
    {
        free(pars->mu_table_inner);
        pars->mu_table_inner = NULL;
    }
}

