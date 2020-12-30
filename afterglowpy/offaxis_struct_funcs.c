#include <string.h>
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

double f_E_exponential2(double theta, void *params)
{
    struct fluxParams *pars = (struct fluxParams *)params;
    if(theta <= pars->theta_wing)
    {
        double x = theta / pars->theta_core;
        double y = theta / 0.225;  // From Lazzati 19
        return pars->E_iso_core * (exp(-x) + pars->b * exp(-y));
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

double check_t_e(double t_e, double mu, double t_obs, double *mu_table, int N)
{
    if(mu > mu_table[N-1])
    {
        fprintf(stderr, "mu >> 1? this should not have happened\n");
        fprintf(stderr,
                "   t_obs=%.6lg t_e=%.6lg mu=%.6lg mu_table[-1]=%.6lg\n",
                t_obs, t_e, mu, mu_table[N-1]);
        return -1;
    }

    if(mu_table[0] >= mu) // happens only if t_e very small
    {
        fprintf(stderr, "very small mu: mu=%.3lg, mu[0]=%.3lg\n",
                mu, mu_table[0]);
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


    double th0 = pars->theta_h;
    double fom = 2*sin(0.5*th0)*sin(0.5*th0); //Fraction of solid angle in jet.

    double Mej_sph;
    if(pars->g_init > 1.0)
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
    double R0, u0;
    //shockInitDecel(Rt0, &R0, &u0, args);
    shockInitFind(Rt0, &R0, &u0, pars->tRes/10, args);
    //printf("t0=%.6le R0=%.6le u0=%.6le\n", Rt0, R0, u0);

    args[0] = pars->E_iso * fom;
    args[1] = Mej_sph * fom;
    shockEvolveSpreadRK4(t_table, R_table, u_table, th_table, table_entries,
                            R0, u0, th0, args, spread);

    if(R_table[0] != R_table[0])
    {
        char msg[MSG_LEN];
        int c = 0;
        c += snprintf(msg, MSG_LEN,
                      "Shock integration Error: R[0]=%.3e  (fac=%.3e)\n",
                      R_table[0], fac);

        c += snprintf(msg+c, MSG_LEN-c,
                      "    t0=%.3e R0=%.3e u0=%.3e th0=%.3e\n",
                      Rt0, R0, u0, th0);
        set_error(pars, msg);
        return;
    }

    if(R_table[table_entries-1] != R_table[table_entries-1])
    {
        char msg[MSG_LEN];
        int c = 0;
        c += snprintf(msg, MSG_LEN,
                      "Shock integration Error: R[-1]=%.3e  (fac=%.3e)\n",
                      R_table[table_entries-1], fac);

        c += snprintf(msg+c, MSG_LEN-c,
                      "    t0=%.3e R0=%.3e u0=%.3e th0=%.3e\n",
                      Rt0, R0, u0, th0);
        set_error(pars, msg);
        return;
    }
}

///////////////////////////////////////////////////////////////////////////////

double emissivity(double nu, double R, double mu, double te,
                    double u, double us, double n0, double p, double epse,
                    double epsB, double ksiN, int specType)
{
    if(us < 1.0e-5)
    {
        //shock is ~ at sound speed of warm ISM. Won't shock, approach invalid.
        return 0.0;
    }
    if(R == 0.0)
        return 0.0;

    // set remaining fluid quantities
    double g = sqrt(1+u*u);
    double beta = u/g;
    double betaS = us / sqrt(1+us*us);
    double nprime = 4.0 * n0 * g; // comoving number density
    double e_th = u*u/(g+1) * nprime * m_p * v_light * v_light;
    double B = sqrt(epsB * 8.0 * PI * e_th);
    double a = (1.0 - mu * beta); // beaming factor
    double ashock = (1.0 - mu * betaS); // shock velocity beaming factor
    double DR = R / (12.0 * g*g * ashock);
    if (DR < 0.0) DR *= -1.0; // DR is function of the absolute value of mu


    double epsebar;
    if(specType & EPS_E_BAR_FLAG)
        epsebar = epse;
    else
        epsebar = (2.0-p) / (1.0-p) * epse;


    // set local emissivity 
    double nuprime = nu * g * a; // comoving observer frequency
    double g_m = epsebar * e_th / (ksiN * nprime * m_e * v_light * v_light);
    double g_c = 6 * PI * m_e * g * v_light / (sigma_T * B * B * te);

    //Inverse Compton adjustment of lfac_c
    if(specType & IC_COOLING_FLAG)
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
    {
        fprintf(stderr, "bad em:%.3le te=%.3le mu=%.3lf\n",
                em, te, mu);
        return -1;
    }
    if(freq != freq || freq < 0.0)
    {
        fprintf(stderr, "bad freq at:%.3le te=%.3le mu=%.3lf\n",
                freq, te, mu);
        return -1;
    }

    double em_lab = em * freq / (g*g * a*a);

    // Self-Absorption
    if(specType & (SSA_SMOOTH_FLAG | SSA_SHARP_FLAG))
    {
        // Co-moving frame absorption coefficient
        double abs_com_P = sqrt(3) * e_e*e_e*e_e * (p-1)*(p+2)*nprime*B
                            / (16*M_PI * m_e*m_e*v_light*v_light
                                * g_m * nuprime*nuprime);
        double abs_com_freq;
        if(nuprime < nu_m)
            abs_com_freq = pow(nuprime / nu_m, 1.0/3.0);
        else
            abs_com_freq = pow(nuprime / nu_m, -0.5*p);

        // Lab frame absorption coefficient
        double abs = abs_com_P * abs_com_freq * a*g;

        // (Signed) Optical depth through this shell.
        // if negative, face is oriented away from observer.
        double dtau;
        if(mu == betaS)
            dtau = 1.0e100; // HUGE VAL, just in case
        else
            dtau = abs * DR * (1 - mu*betaS) / (mu - betaS);


        // Now that we know the optical depth, we apply it in a way
        // according to the given specType

        if((specType & SSA_SMOOTH_FLAG) && (specType & SSA_SHARP_FLAG))
        {
            //Special case: use the optically thick limit *everywhere*
            if(dtau <= 0.0)
                em_lab = 0.0;
            else
                em_lab /= dtau;
        }
        else if(specType & SSA_SMOOTH_FLAG)
        {
            // Apply self-absorption "properly"
            //
            // correction factor to emissivity from absorption
            // ( 1 - e^(-tau) ) / tau  (on front face)
            //
            // back face has extra factor ~ e^-betaS/(mu-betaS)
            //
            // for now ignoring shadowing by the front face.
            double abs_fac;
            if(dtau == 0.0)
                abs_fac = 1.0;
            else if(dtau > 0.0)
                abs_fac = -expm1(-dtau) / dtau;
            else
            {
                abs_fac = expm1(dtau) / dtau; //* exp(
                            //abs * DR * betaS*mu / (mu - betaS));
            }

            em_lab *= abs_fac;
        }
        else if(specType & SSA_SHARP_FLAG)
        {
            // Apply self-absorption "simply".  
            //
            // Compute flux in optically thick limit,
            // use in final result if less than optically thin calculation.
            //
            // e.g. use tau->infty limit if tau > 1.0

            // "Forward" face
            if(dtau > 1.0)
                em_lab /= dtau;
            
            // "Back" face --> assume shadowed by front
            else if(dtau < -1.0)
                em_lab = 0.0;
        }
    }
    if(specType < 0)
        em_lab = 1.0;

    return R * R * DR * em_lab;
}

double costheta_integrand(double aomct, void* params) // inner integral
{
    /*
     * This is the integrand for the inner integral, over theta.
     * The integral is actually performed over 1-cos(theta), 
     * which eliminates the geometrical sin(theta) factor the standard volume
     * element and retains numerical accuracy near theta=0.
     *
     * It is good to know that 1 - cos(theta) = 2*sin(theta/2)^2
     */

    struct fluxParams *pars = (struct fluxParams *) params;

    pars->nevals += 1;
    
    //double cp = cos(pars->phi); 
    //double cto = cos(pars->theta_obs_cur);
    //double sto = sin(pars->theta_obs_cur);
    
    double act = 1 - aomct;
    double a_theta = 2 * asin(sqrt(0.5 * aomct));
    double ast = sqrt(aomct * (1+act));
    pars->theta = a_theta;
    pars->ct = act;
    pars->st = ast;

    //double a_theta = acos(act);
    //double ast = sqrt((1.0 - act)*(1 + act));
    double mu = ast * (pars->cp) * (pars->sto) + act * (pars->cto);

    int ia = searchSorted(mu, pars->mu_table, pars->table_entries);
    int ib = ia+1;
    double t_e = interpolateLin(ia, ib, mu, pars->mu_table, pars->t_table, 
                            pars->table_entries);
    t_e = check_t_e(t_e, mu, pars->t_obs, pars->mu_table, pars->table_entries);

    if(t_e < 0.0)
    {
        char msg[MSG_LEN];
        int c = 0;
        c += snprintf(msg, MSG_LEN-c,
                     "BAD t_e: %.6lf Eiso=%.3le n0=%.3le thetah=%.3le\n",
                     t_e, pars->E_iso, pars->n_0, pars->theta_h);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    theta_obs=%.3lf phi=%.3lf theta=%.3lf mu=%.3lf\n",
                      pars->theta_obs, pars->phi, pars->theta, mu);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    L0=%.3le q=%.3lf ts=%.3le\n",
                      pars->L0, pars->q, pars->ts);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    t[0]=%.3le t[-1]=%.3le R[0]=%.3le R[-1]=%.3le\n",
                      pars->t_table[0], pars->t_table[pars->table_entries-1],
                      pars->R_table[0], pars->R_table[pars->table_entries-1]);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    u[0]=%.3le u[-1]=%.3le th[0]=%.3le th[-1]=%.3le\n",
                      pars->u_table[0], pars->u_table[pars->table_entries-1],
                      pars->th_table[0], pars->th_table[pars->table_entries-1]);
        set_error(pars, msg);
        return 0.0;
    }
    
    double R = interpolateLog(ia, ib, t_e, pars->t_table, pars->R_table, 
                            pars->table_entries);

    double us, u;
    u = interpolateLog(ia, ib, t_e, pars->t_table, pars->u_table,
                                    pars->table_entries);
    us = shockVel(u);
    
    double dFnu =  emissivity(pars->nu_obs, R, mu, t_e, u, us,
                                pars->n_0, pars->p, pars->epsilon_E,
                                pars->epsilon_B, pars->ksi_N, pars->spec_type);

    if(dFnu != dFnu || dFnu < 0.0)
    {
        char msg[MSG_LEN];
        int c = 0;
        c += snprintf(msg, MSG_LEN,
                     "bad dFnu:%.3le nu=%.3le R=%.3le th=%.3lf mu=%.3lf\n",
                     dFnu, pars->nu_obs, R, a_theta, mu);
        c += snprintf(msg+c, MSG_LEN-c,
                     "      t=%.3le u=%.3le us=%.3le n0=%.3le p=%.3lf\n",
                     t_e, u, us, pars->n_0, pars->p);
        c += snprintf(msg+c, MSG_LEN-c,
                     "      epse=%.3le epsB=%.3le ksiN=%.3le specType=%d\n",
                     pars->epsilon_E, pars->epsilon_B, pars->ksi_N,
                     pars->spec_type);
        c += snprintf(msg+c, MSG_LEN-c,
                     "      Rt0=%.3le Rt1=%.3le E_iso=%.3le L0=%.3le ts=%.3le\n",
                     pars->Rt0, pars->Rt1, pars->E_iso, pars->L0, pars->ts);
        set_error(pars, msg);
        return 0.0;
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
    {
        char msg[MSG_LEN];
        snprintf(msg, MSG_LEN, "bad mask fac: %.3le\n", fac);
        set_error(pars, msg);
        return 0.0;
    }

    return fac * dFnu;
}

///////////////////////////////////////////////////////////////////////////////

double phi_integrand(double a_phi, void* params) // outer integral
{
    /*
     * This is the integrand for the (outer) phi integral. This is the
     * inner integral over ~theta.  For stability and smoothness, the
     * integral is performed over 1-cos(theta) instead of over theta 
     * itself.  It is good to know that 1 - cos(theta) = 2 * sin(theta/2)^2
     */

    double result;

    struct fluxParams *pars = (struct fluxParams *) params;
  
    // set up integration routine

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

        if(pars->table_entries_inner == 0)
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
 
    // For a given phi, integrate over 1 - cos(theta)

    double sht0 = sin(0.5*theta_0);
    double sht1 = sin(0.5*theta_1);
    double omct0 = 2 * sht0*sht0;
    double omct1 = 2 * sht1*sht1;

    if(pars->int_type == INT_TRAP_FIXED)
    {
        result = trap(&costheta_integrand, omct0, omct1, pars->nmax_theta,
                      params, check_error);
    }
    else if(pars->int_type == INT_TRAP_ADAPT)
    {
        result = trap_adapt(&costheta_integrand, omct0, omct1,
                            pars->nmax_theta, pars->atol_theta,
                            pars->rtol_theta, params, NULL, NULL, NULL, 0,
                            check_error, NULL, NULL);
    }
    else if(pars->int_type == INT_SIMP_FIXED)
    {
        result = simp(&costheta_integrand, omct0, omct1, pars->nmax_theta,
                      params, check_error);
    }
    else if(pars->int_type == INT_SIMP_ADAPT)
    {
        result = simp_adapt(&costheta_integrand, omct0, omct1,
                            pars->nmax_theta, pars->atol_theta,
                            pars->rtol_theta, params, NULL, NULL, NULL, 0,
                            check_error, NULL, NULL);
    }
    else if(pars->int_type == INT_ROMB_ADAPT)
    {
        int Neval = 0;
        double err = 0;

        result = romb(&costheta_integrand, omct0, omct1, pars->nmax_theta,
                        pars->atol_theta, pars->rtol_theta, params,
                        &Neval, &err, 0, check_error, NULL, NULL);
        //printf("phi = %.3lf:  res=%.6lg  err=%.3lg  Neval=%d  tol=%.3g\n",
        //        a_phi, result, err, Neval,
        //        pars->atol_theta + pars->rtol_theta*result);
    }
    else if(pars->int_type == INT_TRAP_NL)
    {
        result = trapNL_adapt(&costheta_integrand, omct0, omct1,
                              pars->nmax_theta, pars->atol_theta,
                              pars->rtol_theta, params, NULL, NULL, NULL, 0,
                              check_error, NULL, NULL);
    }
    else if(pars->int_type == INT_HYBRID)
    {
        result = hybrid_adapt(&costheta_integrand, omct0, omct1,
                              pars->nmax_theta, pars->atol_theta,
                              pars->rtol_theta, params, NULL, NULL, 0,
                              check_error, NULL, NULL);
    }
    else if(pars->int_type == INT_CADRE)
    {
        result = cadre_adapt(&costheta_integrand, omct0, omct1,
                              pars->nmax_theta, pars->atol_theta,
                              pars->rtol_theta, params, NULL, NULL, 0,
                              check_error, NULL, NULL);
    }
    else if(pars->int_type == INT_GK49_ADAPT)
    {
        result = gk49_adapt(&costheta_integrand, omct0, omct1,
                              pars->nmax_theta, pars->atol_theta,
                              pars->rtol_theta, params, NULL, NULL, 0,
                              check_error);
    }
    else if(pars->int_type == INT_GK715_ADAPT)
    {
        result = gk715_adapt(&costheta_integrand, omct0, omct1,
                              pars->nmax_theta, pars->atol_theta,
                              pars->rtol_theta, params, NULL, NULL, 0,
                              check_error);
    }
    else if(pars->int_type == INT_GK1021_ADAPT)
    {
        result = gk1021_adapt(&costheta_integrand, omct0, omct1,
                              pars->nmax_theta, pars->atol_theta,
                              pars->rtol_theta, params, NULL, NULL, 0,
                              check_error);
    }
    else
    {
        char msg[MSG_LEN];
        snprintf(msg, MSG_LEN,
                 "Unknown integrator %d!  Aborting.\n", pars->int_type);
        set_error(pars, msg);
        return 0.0;
    }
    ERR_CHK_DBL(pars)

    if(result != result || result < 0.0)
    {
        char msg[MSG_LEN];
        int c = 0;
        c += snprintf(msg, MSG_LEN,
                     "bad result in phi_integrand :%.3le\n", result);

        c += snprintf(msg+c, MSG_LEN-c,
                     "   t_obs=%.3le theta_lo=%.3lf theta_hi=%.3lf phi=%.3lf\n",
                     pars->t_obs, theta_0, theta_1, pars->phi);
        set_error(pars, msg);
        return 0.0;
    }
  
    //printf("   a_phi: %.6lf (%.6le)\n", a_phi, result);

    return result;
}

double find_jet_edge(double phi, double cto, double sto, double theta0,
                     double *a_mu, double *a_thj, int N)
{
    /*
     *
     * Return the (outer) edge of the jet section for a particular obs.
     *
     * phi: double
     *      phi-coordinate of jet along which to search
     * cto: double
     *      cos(theta_obs) cosine of observer angle
     * sto: double
     *      sin(theta_obs) sine of observer angle
     * theta0: double
     *      
     * a_mu: double array
     *      time-ordered array of mu values for this observation.
     *      mu = c * (t_em - t_obs) / R(t_em)
     *         = cos(th_obs)*cos(theta) + sin(theta_obs)*sin(theta)*cos(phi)
     * a_thj: double array
     *      time ordered array of jet-edge values.
     */
    
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


///////////////////////////////////////////////////////////////////////////////

double flux(struct fluxParams *pars, double atol) // determine flux for a given t_obs
{
    double result;
    double phi_0 = 0.0;
    double phi_1 = PI;

    // at this stage t_obs is known, so mu_table can be made
    make_mu_table(pars); 

    double d_L = pars->d_L;

    double Fcoeff = cgs2mJy / (4*PI * d_L*d_L);

    //pars->atol_theta = 0.0;
    //double I0 = phi_integrand(0.0, pars);
    //pars->atol_theta = 1.0e-6 * I0;
    
    //Given we want the integral over phi to have an absolute tolerance of
    // atol/(2*Fcoeff), we only need the to know the integrand (the integral
    // over theta) to a tolerance of atol / (2*Fcoeff)
    pars->atol_theta = atol/(2*Fcoeff*PI);

    if(pars->int_type == INT_TRAP_FIXED)
    {
        result = 2 * Fcoeff * trap(&phi_integrand, phi_0, phi_1,
                                    pars->nmax_phi, pars, check_error);
    }
    else if(pars->int_type == INT_TRAP_ADAPT)
    {
        result = 2 * Fcoeff * trap_adapt(&phi_integrand, phi_0, phi_1,
                                         pars->nmax_phi, atol/(2*Fcoeff),
                                         pars->rtol_phi, pars, NULL, NULL,
                                         NULL, 0, check_error, NULL, NULL);
    }
    else if(pars->int_type == INT_SIMP_FIXED)
    {
        result = 2 * Fcoeff * simp(&phi_integrand, phi_0, phi_1,
                                    pars->nmax_phi, pars, check_error);
    }
    else if(pars->int_type == INT_SIMP_ADAPT)
    {
        result = 2 * Fcoeff * simp_adapt(&phi_integrand, phi_0, phi_1,
                                         pars->nmax_phi, atol/(2*Fcoeff),
                                         pars->rtol_phi, pars, NULL, NULL,
                                         NULL, 0, check_error, NULL, NULL);
    }
    else if(pars->int_type == INT_ROMB_ADAPT)
    {
        double phi_a = phi_0 + 0.5*(phi_1-phi_0);
        result = 2 * Fcoeff * romb(&phi_integrand, phi_0, phi_a,
                                    pars->nmax_phi, atol/(2*Fcoeff),
                                    pars->rtol_phi, pars, NULL, NULL, 0,
                                    check_error, NULL, NULL);
        ERR_CHK_DBL(pars)
        result += 2 * Fcoeff * romb(&phi_integrand, phi_a, phi_1,
                                    pars->nmax_phi,
                                    (atol+pars->rtol_phi*result)/(2*Fcoeff),
                                    pars->rtol_phi, pars, NULL, NULL, 0,
                                    check_error, NULL, NULL);
    }
    else if(pars->int_type == INT_TRAP_NL)
    {
        result = 2 * Fcoeff * trapNL_adapt(&phi_integrand, phi_0, phi_1,
                                           pars->nmax_phi, atol/(2*Fcoeff),
                                           pars->rtol_phi, pars, NULL, NULL,
                                           NULL, 0, check_error, NULL, NULL);
    }
    else if(pars->int_type == INT_HYBRID)
    {
        result = 2 * Fcoeff * hybrid_adapt(&phi_integrand, phi_0, phi_1,
                                           pars->nmax_phi, atol/(2*Fcoeff),
                                           pars->rtol_phi, pars, NULL, NULL,
                                           0, check_error, NULL, NULL);
    }
    else if(pars->int_type == INT_CADRE)
    {
        result = 2 * Fcoeff * cadre_adapt(&phi_integrand, phi_0, phi_1,
                                           pars->nmax_phi, atol/(2*Fcoeff),
                                           pars->rtol_phi, pars, NULL, NULL,
                                           0, check_error, NULL, NULL);
    }
    else if(pars->int_type == INT_GK49_ADAPT)
    {
        result = 2 * Fcoeff * gk49_adapt(&phi_integrand, phi_0, phi_1,
                                           pars->nmax_phi, atol/(2*Fcoeff),
                                           pars->rtol_phi, pars, NULL, NULL,
                                           0, check_error);
    }
    else if(pars->int_type == INT_GK715_ADAPT)
    {
        result = 2 * Fcoeff * gk715_adapt(&phi_integrand, phi_0, phi_1,
                                           pars->nmax_phi, atol/(2*Fcoeff),
                                           pars->rtol_phi, pars, NULL, NULL,
                                           0, check_error);
    }
    else if(pars->int_type == INT_GK1021_ADAPT)
    {
        result = 2 * Fcoeff * gk1021_adapt(&phi_integrand, phi_0, phi_1,
                                           pars->nmax_phi, atol/(2*Fcoeff),
                                           pars->rtol_phi, pars, NULL, NULL,
                                           0, check_error);
    }
    else
    {
        char msg[MSG_LEN];
        snprintf(msg, MSG_LEN,
                 "Unknown integrator %d!  Aborting.\n", pars->int_type);
        set_error(pars, msg);
        return 0.0;
    }

    ERR_CHK_DBL(pars)
    
    if(result != result || result < 0.0)
    {
        char msg[MSG_LEN];
        int c = 0;
        c += snprintf(msg, MSG_LEN,
                     "bad result in flux() :%.3le\n", result);

        c += snprintf(msg+c, MSG_LEN-c,
                "   t_obs=%.3le nu_obs=%.3le theta_lo=%.3lf theta_hi=%.3lf\n",
                     pars->t_obs, pars->nu_obs, pars->current_theta_cone_low,
                     pars->current_theta_cone_hi);
        c += snprintf(msg+c, MSG_LEN-c,
                "   Fcoeff=%.6le\n", Fcoeff);
        set_error(pars, msg);
        return 0.0;
    }

    return result;
}

void lc_cone(double *t, double *nu, double *F, int Nt, double E_iso,
                double theta_core, double theta_wing, struct fluxParams *pars)
{
    int i;

    set_jet_params(pars, E_iso, theta_wing);
    ERR_CHK_VOID(pars)

    for(i=0; i<Nt; i++)
    {
        F[i] = flux_cone(t[i], nu[i], -1, -1, theta_core, theta_wing, 0.0,
                            pars);
        ERR_CHK_VOID(pars)
    }
}

void lc_tophat(double *t, double *nu, double *F, int Nt,
                double E_iso, double theta_h, struct fluxParams *pars)
{
    int i;

    set_jet_params(pars, E_iso, theta_h);
    ERR_CHK_VOID(pars)

    for(i=0; i<Nt; i++)
    {
        F[i] = flux_cone(t[i], nu[i], -1, -1, 0.0, theta_h, 0.0, pars);
        ERR_CHK_VOID(pars)
    }
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

        //printf("cone %d: th_lo=%.6lf th_hi=%.6lf, E=%.6le\n", i,
        //        theta_cone_low, theta_cone_hi, E_iso);

        if(theta_c_arr != NULL)
            theta_c_arr[i] = theta_c;
        if(E_iso_arr != NULL)
            E_iso_arr[i] = E_iso;

        if(E_iso <= 0.0)
            continue;

        set_jet_params(pars, E_iso, theta_h);
        ERR_CHK_VOID(pars)

        for(j=0; j<Nt; j++)
        {
            //printf("tobs = %.6le\n", t[j]);
            F[j] += flux_cone(t[j], nu[j], -1, -1, theta_cone_low,
                                theta_cone_hi,
                                F[j]*pars->rtol_struct/res_cones,
                                pars);
            ERR_CHK_VOID(pars)
        }
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
    ERR_CHK_VOID(pars)

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

        if(E_iso <= 0.0)
            continue;

        set_jet_params(pars, E_iso, theta_h);
        ERR_CHK_VOID(pars)

        for(j=0; j<Nt; j++)
        {
            F[j] += flux_cone(t[j], nu[j], -1, -1, theta_cone_low,
                                theta_cone_hi,
                                F[j]*pars->rtol_struct/res_cones,
                                pars);
            ERR_CHK_VOID(pars)
        }
    }
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
    {
        set_jet_params(pars, E_iso, theta_h);
        ERR_CHK_DBL(pars)
    }

    //Jet 
    theta_obs_cur = theta_obs;
    set_obs_params(pars, t_obs, nu_obs, theta_obs_cur, 
                    theta_cone_hi, theta_cone_low);
    F1 = flux(pars, atol);
    ERR_CHK_DBL(pars)
    
    //Counter-jet
    if(pars->counterjet)
    {
        theta_obs_cur = 180*deg2rad - theta_obs;
        set_obs_params(pars, t_obs, nu_obs, theta_obs_cur, 
                        theta_cone_hi, theta_cone_low);
        F2 = flux(pars, atol);
        ERR_CHK_DBL(pars)
    }
    else
        F2 = 0.0;
    
    Fboth = F1 + F2;


    if(F1 != F1 || F1 < 0.0)
    {
        char msg[MSG_LEN];
        int c = snprintf(msg, MSG_LEN, "bad F1 in flux_cone:%.3lg\n", F1);
        c += snprintf(msg+c, MSG_LEN-c,
                      "      t_obs=%.3le theta_lo=%.3lf theta_hi=%.3lf\n",
                      t_obs, theta_cone_low, theta_cone_hi);
        set_error(pars, msg);
        return 0.0;
    }
    if(F2 != F2 || F2 < 0.0)
    {
        char msg[MSG_LEN];
        int c = snprintf(msg, MSG_LEN, "bad F2 in flux_cone:%.3lg\n", F2);
        c += snprintf(msg+c, MSG_LEN-c,
                      "      t_obs=%.3le theta_lo=%.3lf theta_hi=%.3lf\n",
                      t_obs, theta_cone_low, theta_cone_hi);
        set_error(pars, msg);
        return 0.0;
    }

    //printf(" Fcone = %.6le\n", Fboth);

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
    {
        char msg[MSG_LEN];
        int c = 0;
        c += snprintf(msg, MSG_LEN-c,
                     "BAD t_e: %.6lf Eiso=%.3le n0=%.3le thetah=%.3le\n",
                     t_e, pars->E_iso, pars->n_0, pars->theta_h);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    theta_obs=%.3lf phi=%.3lf theta=%.3lf mu=%.3lf\n",
                      pars->theta_obs, pars->phi, pars->theta, mu);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    L0=%.3le q=%.3lf ts=%.3le\n",
                      pars->L0, pars->q, pars->ts);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    t[0]=%.3le t[-1]=%.3le R[0]=%.3le R[-1]=%.3le\n",
                      pars->t_table[0], pars->t_table[pars->table_entries-1],
                      pars->R_table[0], pars->R_table[pars->table_entries-1]);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    u[0]=%.3le u[-1]=%.3le th[0]=%.3le th[-1]=%.3le\n",
                      pars->u_table[0], pars->u_table[pars->table_entries-1],
                      pars->th_table[0], pars->th_table[pars->table_entries-1]);
        set_error(pars, msg);
        return 0.0;
    }

    double R = interpolateLog(ia, ib, t_e, pars->t_table,
                                pars->R_table, pars->table_entries);
    double u = interpolateLog(ia, ib, t_e, pars->t_table,
                                pars->u_table, pars->table_entries);
    double us = shockVel(u);

    I = emissivity(pars->nu_obs, R, mu, t_e, u, us, pars->n_0,
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
    {
        char msg[MSG_LEN];
        int c = 0;
        c += snprintf(msg, MSG_LEN-c,
                     "BAD t_e: %.6lf Eiso=%.3le n0=%.3le thetah=%.3le\n",
                     t_e, pars->E_iso, pars->n_0, pars->theta_h);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    theta_obs=%.3lf phi=%.3lf theta=%.3lf mu=%.3lf\n",
                      pars->theta_obs, pars->phi, pars->theta, mu);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    L0=%.3le q=%.3lf ts=%.3le\n",
                      pars->L0, pars->q, pars->ts);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    t[0]=%.3le t[-1]=%.3le R[0]=%.3le R[-1]=%.3le\n",
                      pars->t_table[0], pars->t_table[pars->table_entries-1],
                      pars->R_table[0], pars->R_table[pars->table_entries-1]);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    u[0]=%.3le u[-1]=%.3le th[0]=%.3le th[-1]=%.3le\n",
                      pars->u_table[0], pars->u_table[pars->table_entries-1],
                      pars->th_table[0], pars->th_table[pars->table_entries-1]);
        set_error(pars, msg);
        return;
    }

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
    ERR_CHK_VOID(pars)
    set_obs_params(pars, t[0], nu[0], theta_obs,
                   theta_cone_hi, theta_cone_low);
    make_mu_table(pars);
    double tobs = t[0];

    for(j=0; j<N; j++)
    {
        double th = theta[j];
        double ph = phi[j];

        if(I[j] > 0.0 || th < theta_cone_low)
            continue;

        set_obs_params(pars, t[j], nu[j], theta_obs,
                       theta_cone_hi, theta_cone_low);
        if(tobs != t[j])
        {
            make_mu_table(pars);
            tobs = t[j];
        }

        double th_a, th_b;
        th_b = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                             pars->mu_table, pars->th_table,
                             pars->table_entries);
        if(pars->table_entries_inner == 0)
            th_a = (theta_cone_low / theta_cone_hi) * th_b;
        else
            th_a = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_low,
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
        ERR_CHK_VOID(pars)

        set_obs_params(pars, t[0], nu[0], theta_obs, theta_cone_hi, 
                       theta_cone_low);
        make_mu_table(pars);
        double tobs = t[0];

        //printf("Cone: %d of %d\n", i, res_cones);

        for(j=0; j<N; j++)
        {
            double th = theta[j];
            double ph = phi[j];

            if(I[j] > 0.0 || th < theta_cone_low)
                continue;


            set_obs_params(pars, t[j], nu[j], theta_obs, theta_cone_hi, 
                            theta_cone_low);
            if(tobs != t[j])
            {
                make_mu_table(pars);
                tobs = t[j];
            }

            double th_a, th_b;
            th_b = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                 pars->mu_table, pars->th_table,
                                 pars->table_entries);
            ERR_CHK_VOID(pars)
            if(pars->table_entries_inner == 0)
                th_a = (theta_cone_low / theta_cone_hi) * th_b;
            else
            {
                th_a = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_low,
                                     pars->mu_table_inner, 
                                     pars->th_table_inner,
                                     pars->table_entries_inner);
                ERR_CHK_VOID(pars)
            }

            if(th < th_a || th > th_b)
                continue;

            I[j] += Fcoeff * intensity(th, ph, t[j], nu[j], theta_obs,        
                                    theta_cone_hi, theta_cone_low, pars);
            ERR_CHK_VOID(pars)
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
    ERR_CHK_VOID(pars)

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
        ERR_CHK_VOID(pars)
        set_obs_params(pars, t[0], nu[0], theta_obs, theta_cone_hi, 
                        theta_cone_low);
        make_mu_table(pars);
        double tobs = t[0];
        
        for(j=0; j<N; j++)
        {
            double th = theta[j];
            double ph = phi[j];
            
            if(I[j] > 0.0 || th < theta_cone_low)
                continue;

            set_obs_params(pars, t[j], nu[j], theta_obs, theta_cone_hi, 
                            theta_cone_low);
            if(tobs != t[j])
            {
                make_mu_table(pars);
                tobs = t[j];
            }

            double th_a, th_b;
            th_b = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                 pars->mu_table, pars->th_table,
                                 pars->table_entries);
            ERR_CHK_VOID(pars)
            if(pars->table_entries_inner == 0)
                th_a = (theta_cone_low / theta_cone_hi) * th_b;
            else
            {
                th_a = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_low,
                                     pars->mu_table_inner, 
                                     pars->th_table_inner,
                                     pars->table_entries_inner);
                ERR_CHK_VOID(pars)
            }

            if(th < th_a || th > th_b)
                continue;

            I[j] += Fcoeff * intensity(th, ph, t[j], nu[j], theta_obs,        
                                    theta_cone_hi, theta_cone_low, pars);
            ERR_CHK_VOID(pars)
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
    ERR_CHK_VOID(pars)
    set_obs_params(pars, tobs[0], 1.0, theta_obs, theta_cone_hi, 
                    theta_cone_low);
    make_mu_table(pars);
    double tobs_cur = t[0];

    for(j=0; j<N; j++)
    {
        double th = theta[j];
        double ph = phi[j];

        if(t[j] > 0.0 || th < theta_cone_low)
            continue;

        set_obs_params(pars, tobs[j], 1.0, theta_obs,
                       theta_cone_hi, theta_cone_low);
        if(tobs_cur != tobs[j])
        {
            make_mu_table(pars);
            tobs_cur = tobs[j];
        }

        double th_a, th_b;
        th_b = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                             pars->mu_table, pars->th_table,
                             pars->table_entries);
        ERR_CHK_VOID(pars)
        if(pars->table_entries_inner == 0)
            th_a = (theta_cone_low / theta_cone_hi) * th_b;
        else
        {
            th_a = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_low,
                                 pars->mu_table_inner, 
                                 pars->th_table_inner,
                                 pars->table_entries_inner);
            ERR_CHK_VOID(pars)
        }

        if(th < th_a || th > th_b)
            continue;

        shockVals(th, ph, tobs[j], t+j, R+j, u+j, thj+j,
                    theta_obs, theta_cone_hi, theta_cone_low, pars);
        ERR_CHK_VOID(pars)
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
        ERR_CHK_VOID(pars)
        set_obs_params(pars, tobs[0], 1.0, theta_obs, theta_cone_hi, 
                        theta_cone_low);
        make_mu_table(pars);
        double tobs_cur = t[0];

        for(j=0; j<N; j++)
        {
            double th = theta[j];
            double ph = phi[j];

            if(t[j] > 0.0 || th < theta_cone_low)
                continue;

            set_obs_params(pars, tobs[j], 1.0, theta_obs, theta_cone_hi, 
                            theta_cone_low);
            if(tobs_cur != tobs[j])
            {
                make_mu_table(pars);
                tobs_cur = tobs[j];
            }

            double th_a, th_b;
            th_b = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                 pars->mu_table, pars->th_table,
                                 pars->table_entries);
            ERR_CHK_VOID(pars)
            if(pars->table_entries_inner == 0)
                th_a = (theta_cone_low / theta_cone_hi) * th_b;
            else
            {
                th_a = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_low,
                                     pars->mu_table_inner, 
                                     pars->th_table_inner,
                                     pars->table_entries_inner);
                ERR_CHK_VOID(pars)
            }

            if(th < th_a || th > th_b)
                continue;

            shockVals(th, ph, tobs[j], t+j, R+j, u+j, thj+j,
                      theta_obs, theta_cone_hi, theta_cone_low, pars);
            ERR_CHK_VOID(pars)
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
    ERR_CHK_VOID(pars)
            

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
        ERR_CHK_VOID(pars)
        set_obs_params(pars, tobs[0], 1.0, theta_obs, theta_cone_hi, 
                        theta_cone_low);
        make_mu_table(pars);
        double tobs_cur = tobs[0];
        
        for(j=0; j<N; j++)
        {
            double th = theta[j];
            double ph = phi[j];
            
            if(t[j] > 0.0 || th < theta_cone_low)
                continue;

            set_obs_params(pars, tobs[j], 1.0, theta_obs, theta_cone_hi, 
                            theta_cone_low);
            if(tobs_cur != tobs[j])
            {
                make_mu_table(pars);
                tobs_cur = tobs[j];
            }

            double th_a, th_b;
            th_b = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_hi,
                                 pars->mu_table, pars->th_table,
                                 pars->table_entries);
            ERR_CHK_VOID(pars)
            if(pars->table_entries_inner == 0)
                th_a = (theta_cone_low / theta_cone_hi) * th_b;
            else
            {
                th_a = find_jet_edge(ph, pars->cto, pars->sto, theta_cone_low,
                                     pars->mu_table_inner, 
                                     pars->th_table_inner,
                                     pars->table_entries_inner);
                ERR_CHK_VOID(pars)
            }

            if(th < th_a || th > th_b)
                continue;

            shockVals(th, ph, tobs[j], t+j, R+j, u+j, thj+j,
                      theta_obs, theta_cone_hi, theta_cone_low, pars);
            ERR_CHK_VOID(pars)
        }
    }
}

void calc_flux_density(int jet_type, int spec_type, double *t, double *nu,
                            double *Fnu, int N, struct fluxParams *fp)
{
    int latRes = fp->latRes;
    double E_iso_core = fp->E_iso_core;
    double theta_h_core = fp->theta_core;
    double theta_h_wing = fp->theta_wing;

    int res_cones = (int) (latRes*theta_h_wing / theta_h_core);


    if(jet_type == _tophat)
    {
        lc_tophat(t, nu, Fnu, N, E_iso_core, theta_h_core, fp);
    }
    else if(jet_type == _cone)
    {
        lc_cone(t, nu, Fnu, N, E_iso_core, theta_h_core, theta_h_wing, fp);
    }
    else if(jet_type == _Gaussian)
    {
        lc_struct(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_Gaussian, fp);
    }
    else if(jet_type == _powerlaw)
    {
        lc_struct(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_powerlaw, fp);
    }
    else if(jet_type == _twocomponent)
    {
        lc_struct(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_twocomponent, fp);
    }
    else if(jet_type == _Gaussian_core)
    {
        lc_structCore(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_GaussianCore, fp);
    }
    else if(jet_type == _powerlaw_core)
    {
        lc_structCore(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_powerlawCore, fp);
    }
    else if(jet_type == _exponential)
    {
        lc_structCore(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_exponential, fp);
    }
    else if(jet_type == _exponential2)
    {
        lc_structCore(t, nu, Fnu, N, E_iso_core, theta_h_core, 
                theta_h_wing, NULL, NULL, res_cones, &f_E_exponential2, fp);
    }

    //printf("  Calc took %ld evalutions\n", fp->nevals);
}    

void calc_intensity(int jet_type, int spec_type, double *theta, double *phi,
                            double *t, double *nu, double *Inu, int N,
                            struct fluxParams *fp)
{
    int latRes = fp->latRes;
    double E_iso_core = fp->E_iso_core;
    double theta_h_core = fp->theta_core;
    double theta_h_wing = fp->theta_wing;

    int res_cones = (int) (latRes*theta_h_wing / theta_h_core);


    if(jet_type == _tophat)
    {
        intensity_cone(theta, phi, t, nu, Inu, N, E_iso_core, 0.0, theta_h_core,
                            fp);
    }
    else if(jet_type == _cone)
    {
        intensity_cone(theta, phi, t, nu, Inu, N, E_iso_core, theta_h_core,
                        theta_h_wing, fp);
    }
    else if(jet_type == _Gaussian)
    {
        intensity_struct(theta, phi, t, nu, Inu, N, E_iso_core, theta_h_core, 
                theta_h_wing, res_cones, &f_E_Gaussian, fp);
    }
    else if(jet_type == _powerlaw)
    {
        intensity_struct(theta, phi, t, nu, Inu, N, E_iso_core, theta_h_core, 
                theta_h_wing, res_cones, &f_E_powerlaw, fp);
    }
    else if(jet_type == _Gaussian_core)
    {
        intensity_structCore(theta, phi, t, nu, Inu, N, E_iso_core, 
                            theta_h_core, theta_h_wing, res_cones,
                            &f_E_GaussianCore, fp);
    }
    else if(jet_type == _powerlaw_core)
    {
        intensity_structCore(theta, phi, t, nu, Inu, N, E_iso_core, 
                                theta_h_core, theta_h_wing, res_cones, 
                                &f_E_powerlawCore, fp);
    }

}

void calc_shockVals(int jet_type, double *theta, double *phi, double *tobs,
                    double *t, double *R, double *u, double *thj, int N,
                    struct fluxParams *fp)
{
    int latRes = fp->latRes;
    double E_iso_core = fp->E_iso_core;
    double theta_h_core = fp->theta_core;
    double theta_h_wing = fp->theta_wing;

    int res_cones = (int) (latRes*theta_h_wing / theta_h_core);


    if(jet_type == _tophat)
    {
        shockVals_cone(theta, phi, tobs, t, R, u, thj, N, E_iso_core, 0.0,
                        theta_h_core, fp);
    }
    else if(jet_type == _cone)
    {
        shockVals_cone(theta, phi, tobs, t, R, u, thj, N, E_iso_core,
                        theta_h_core, theta_h_wing, fp);
    }
    else if(jet_type == _Gaussian)
    {
        shockVals_struct(theta, phi, tobs, t, R, u, thj, N, E_iso_core,
                         theta_h_core, theta_h_wing, res_cones, &f_E_Gaussian,
                         fp);
    }
    else if(jet_type == _powerlaw)
    {
        shockVals_struct(theta, phi, tobs, t, R, u, thj, N, E_iso_core,
                         theta_h_core, theta_h_wing, res_cones, &f_E_powerlaw,
                         fp);
    }
    else if(jet_type == _Gaussian_core)
    {
        shockVals_structCore(theta, phi, tobs, t, R, u, thj, N, E_iso_core,
                             theta_h_core, theta_h_wing, res_cones,
                             &f_E_GaussianCore, fp);
    }
    else if(jet_type == _powerlaw_core)
    {
        shockVals_structCore(theta, phi, tobs, t, R, u, thj, N, E_iso_core,
                             theta_h_core, theta_h_wing, res_cones,
                             &f_E_powerlawCore, fp);
    }

}
///////////////////////////////////////////////////////////////////////////////

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
                    int spread, int counterjet, int gamma_type)
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
    pars->latRes = latRes;
    pars->int_type = int_type;
    pars->rtol_struct = rtol_struct;
    pars->rtol_phi = rtol_phi;
    pars->rtol_theta = rtol_theta;
    pars->nmax_phi = nmax_phi;
    pars->nmax_theta = nmax_theta;
    
    pars->atol_theta = 0.0;

    pars->mask = mask;
    pars->nmask = nmask;
    pars->spread = spread;
    pars->counterjet = counterjet;

    pars->nevals = 0;

    pars->error = 0;
    pars->error_msg = NULL;
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

    if(pars->gamma_type == GAMMA_FLAT)
        pars->g_init = pars->g_core;
    else if(pars->gamma_type == GAMMA_EVENMASS)
    {
        if(pars->E_core_global > 0.0)
            pars->g_init = 1.0 + (pars->g_core - 1)
                                    * E_iso / pars->E_core_global;
        else
            pars->g_init = 1.0 + (pars->g_core - 1) * E_iso / pars->E_iso_core;
    }
    else
        pars->g_init = -1;

    if(pars->g_init <= 1.0 && (pars->gamma_type == GAMMA_FLAT
        || pars->gamma_type == GAMMA_EVENMASS
        || pars->gamma_type == GAMMA_STRUCT))
    {
        char msg[MSG_LEN];
        int c = 0;
        c += snprintf(msg, MSG_LEN,
                      "Bad initial Lorentz factor: gamma_init = %.6lg\n",
                      pars->g_init);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    E_iso = %.6lg   theta_h = %.6lg\n",
                      E_iso, theta_h);
        c += snprintf(msg+c, MSG_LEN-c,
                      "    gamma_type = %d   gamma_core = %.6lg\n",
                      pars->gamma_type, pars->g_core);
        set_error(pars, msg);
        return;
    }


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

int check_error(void *params)
{
    struct fluxParams *fp = (struct fluxParams *) params;
    return fp->error;
}

///////////////////////////////////////////////////////////////////////////////

void set_error(struct fluxParams *pars, char msg[])
{
    fprintf(stderr, "error in afterglowpy\n");
    pars->error = 1;

    int msglen = (int)strlen(msg);
    char dump[DUMP_MSG_LEN_MAX];
    int c = snprintf(dump, DUMP_MSG_LEN_MAX, "fluxParamsDump\n{\n");
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    theta: %.12lg\n",
            pars->theta);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    phi: %.12lg\n", pars->phi);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    cp: %.12lg\n", pars->cp);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    ct: %.12lg\n", pars->ct);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    st: %.12lg\n", pars->st);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    cto: %.12lg\n", pars->cto);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    sto: %.12lg\n", pars->sto);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    theta_obs: %.12lg\n",
            pars->theta_obs);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    t_obs: %.12lg\n",
            pars->t_obs);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    nu_obs: %.12lg\n",
            pars->nu_obs);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    d_L: %.12lg\n", pars->d_L);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    E_iso: %.12lg\n",
            pars->E_iso);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    n_0: %.12lg\n", pars->n_0);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    g_init: %.12lg\n",
            pars->g_init);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    p: %.12lg\n", pars->p);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    epsilon_E: %.12lg\n",
                pars->epsilon_E);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    epsilon_B: %.12lg\n",
            pars->epsilon_B);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    ksi_N: %.12lg\n",
            pars->ksi_N);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    theta_h: %.12lg\n",
            pars->theta_h);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    E_iso_core: %.12lg\n",
            pars->E_iso_core);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    theta_core: %.12lg\n",
            pars->theta_core);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    theta_wing: %.12lg\n",
            pars->theta_wing);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    b: %.12lg\n", pars->b);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    E_tot: %.12lg\n",
            pars->E_tot);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    g_core: %.12lg\n",
            pars->g_core);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    E_core_global: %.12lg\n",
            pars->E_core_global);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c,
            "    theta_core_global: %.12lg\n", pars->theta_core_global);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    envType: %d\n",
            pars->envType);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    As: %.12lg\n", pars->As);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    Rwind: %.12lg\n",
            pars->Rwind);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    L0: %.12lg\n", pars->L0);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    q: %.12lg\n", pars->q);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    ts: %.12lg\n", pars->ts);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c,
            "    current_theta_cone_hi: %.12lg\n",
            pars->current_theta_cone_hi);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c,
            "    current_theta_cone_low: %.12lg\n", 
            pars->current_theta_cone_low);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    theta_obs_cur: %.12lg\n",
            pars->theta_obs_cur);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    tRes: %d\n", pars->tRes);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    latRes: %d\n",
            pars->latRes);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    spread: %d\n",
            pars->spread);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    counterjet: %d\n",
            pars->counterjet);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    int_type: %d\n",
            pars->int_type);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    rtol_struct: %.12lg\n",
            pars->rtol_struct);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    rtol_theta: %.12lg\n",
            pars->rtol_theta);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    rtol_phi: %.12lg\n",
            pars->rtol_phi);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    nmax_theta: %d\n",
            pars->nmax_theta);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    nmax_phi: %d\n",
            pars->nmax_phi);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    atol_theta: %.12lg\n",
            pars->atol_theta);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    Rt0: %.12lg\n", pars->Rt0);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    Rt1: %.12lg\n", pars->Rt1);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    ta: %.12lg\n", pars->ta);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    tb: %.12lg\n", pars->tb);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    C_BMsqrd: %.12lg\n",
            pars->C_BMsqrd);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    C_STsqrd: %.12lg\n",
            pars->C_STsqrd);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    t_NR: %.12lg\n",
            pars->t_NR);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    table_entries: %d\n",
            pars->table_entries);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    table_entries_inner: %d\n",
            pars->table_entries_inner);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    spec_type: %d\n",
            pars->spec_type);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    gamma_type: %d\n",
            pars->gamma_type);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    nmask: %d\n",
            pars->nmask);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    nevals: %ld\n",
            pars->nevals);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "    error: %d\n", pars->error);
    c += snprintf(dump+c, DUMP_MSG_LEN_MAX-c, "}\n");

    int dumplen = strlen(dump);
    int len = msglen + dumplen + 1;

    pars->error_msg = (char *)malloc(len * sizeof(char));
    
    strcpy(pars->error_msg, msg);
    strcpy(pars->error_msg+msglen, dump);
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

    if(pars->error_msg != NULL)
    {
        free(pars->error_msg);
        pars->error_msg = NULL;
    }
}

