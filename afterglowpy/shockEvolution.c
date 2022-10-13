#include <stdio.h>
#include <math.h>
#include "offaxis_struct.h"
#include "shockEvolution.h"


double shockVel(double u)
{
    return 4*u*sqrt((u*u+1) / (8*u*u+9));
}

double E_inj(double te, double L0, double q, double ts)
{
    if(te < T0_inj)
        return L0*te;

    if (te > ts)
        te = ts;
    
    double E = L0*T0_inj;

    if(q == 0.0)
        return L0*te;
    else if(q == 1.0)
        return E + L0*T0_inj*log(te/T0_inj);
    else
        return E + L0*T0_inj*(pow(te/T0_inj,1-q) - 1.0)/(1-q);

    return 0.0;
}

double L_inj(double te, double L0, double q, double ts)
{
    if(te <= T0_inj)
        return L0;
    else if(te < ts)
    {
        if(q == 0.0)
            return L0;
        else
            return L0*pow(te/T0_inj,-q);
    }

    return 0.0;
}

double te_inj(double Ei, double L0, double q, double ts)
{
    if(L0*T0_inj >= Ei)
        return Ei / L0;
    Ei -= L0*T0_inj;

    double te;

    if(q == 0)
        te = Ei/L0;
    else if (q == 1)
        te = T0_inj * exp(Ei/(L0*T0_inj));
    else
        te = T0_inj * pow((1-q)*Ei/(L0*T0_inj)+1, 1.0/(1-q)); 

    if(te > ts)
        return -1.0;

    return te;
}

double t_inj(double te, double t0, double te0, double L0, double q, double ts)
{
    return t0*pow(1 + 16*(te-te0)/((2-q)*t0), 0.25*(2-q));
}


double envDensityPar(double R, struct shockParams *par)
{
    return envDensity(R, par->envType, par->rho0_env, par->R0_env, par->k_env,
                      par->rho1_env);
}

double envMassPar(double R, struct shockParams *par)
{
    return envMass(R, par->envType, par->rho0_env, par->R0_env, par->k_env,
                   par->rho1_env);
}

double envRadiusPar(double M, struct shockParams *par)
{
    return envRadius(M, par->envType, par->rho0_env, par->R0_env, par->k_env,
                     par->rho1_env);
}

double envDensity(double R, int envType, double rho0, double R0, double k,
                  double rho1)
{
    if(envType == ENV_ISM)
        return rho0;
    else if(envType == ENV_WIND)
        return rho0 * R0*R0/(R*R);
    else if(envType == ENV_PL)
        return rho0 * pow(R/R0, -k);
    else if(envType == ENV_STEP)
        return R < R0 ? rho0 : rho1;
    else
        return 0;
}

double envMass(double R, int envType, double rho0, double R0, double k,
               double rho1)
{
    if(envType == ENV_ISM)
        return 4.0/3.0 * M_PI * rho0 * R*R*R;
    else if(envType == ENV_WIND)
        return 4.0*M_PI * rho0 * R0*R0 * R;
    else if(envType == ENV_PL)
        return 4.0*M_PI / (3.0-k) * rho0 * pow(R/R0, 3-k) * R0*R0*R0;
    else if (envType == ENV_STEP)
    {
        if(R < R0)
            return 4.0/3.0 * M_PI * rho0 * R*R*R;
        else
            return 4.0/3.0 * M_PI * (rho0 * R0*R0*R0
                                     + rho1 * (R-R0)*(R-R0)*(R-R0));
    }
    else
        return 0;
}

double envRadius(double M, int envType, double rho0, double R0, double k,
                 double rho1)
{
    if(envType == ENV_ISM)
        return pow(0.75 * M / (M_PI * rho0), 1.0/3.0);
    else if(envType == ENV_WIND)
        return M / (4.0*M_PI * rho0 * R0*R0);
    else if(envType == ENV_PL)
        return R0 * pow((3.0-k) * M / (4.0*M_PI * rho0 * R0*R0*R0), 1.0/(3-k));
    else if(envType == ENV_STEP)
    {
        double M0 = 4.0/3.0 * M_PI * rho0 * R0*R0*R0;
        if(M <= M0)
            return pow(0.75 * M / (M_PI * rho0), 1.0/3.0);
        else
            return R0 + pow(0.75 * (M-M0) / (M_PI * rho1), 1.0/3.0);
    }
    else
        return 0;
}


void shockInitDecel(double t0, double *R0, double *u0, void *argv)
{
    double *args = (double *)argv;
    double E0 = args[0];
    double rho0 = args[2];
    double L0 = args[6];
    double q = args[7];
    double ts = args[8];

    double dr, R, u;
    double c = v_light;
    double c5 = c*c*c*c*c;

    //First pass: assume no energy injection
    double C = sqrt(9/(16.0*M_PI) * E0/(rho0*c5));
    u = C * pow(t0, -1.5);
    dr = 1.0/(u*u*16.0);
    R = c*t0*(1-dr);

    if(L0 < 0.0 || ts < 0.0)
    {
        //printf("No energy injection! E=%.3le\n", E0);
        *R0 = R;
        *u0 = u;
        return;
    }

    //Check if energy injection is important
    double te = t0*dr;
    double Ei = E_inj(te, L0, q, ts);

    if(Ei <= E0)
    {
        //printf("Energy injection not important! E0=%.3le Ei=%.3le\n", E0, Ei);
        *R0 = R;
        *u0 = u;
        return;
    }

    //Second pass: use energy injection solution.

    //Time for transition
    double teb = te_inj(E0, L0, q, ts);
    double t1 = pow(16*C*C*teb, 0.25);

    //position at transition
    double u1 = C * pow(t1, -1.5);
    double dr1 = 1.0/(u1*u1*16.0);

    //Account for initial constant Luminosity
    if(teb < T0_inj)
    {
        //Lab time for T0_inj
        double t0i = t_inj(T0_inj, t1, teb, L0, 0.0, ts);

        if(t0 < t0i) //Praise be
        {
            //printf("Early energy injection.\n");
            u = u1*pow(t0/t1, -0.5);
            dr = (t1*dr1 + 2*(t0/(u*u)-t1/(u1*u1))/16.0) / t0;
            *R0 = c*t0*(1-dr);
            *u0 = u;
            return;
        }
        
        //Move solution forwards.
        u = u1*pow(t0i/t1, -0.5);
        dr = (t1*dr1 + 2*(t0i/(u*u)-t1/(u1*u1))/16.0) / t0i;

        t1 = t0i;
        u1 = u;
        dr1 = dr;
    }

    double te1 = dr1*t1;

    //Get lab time when energy injection ends.
    double tls = t_inj(ts, t1, te1, L0, q, ts);
    if(t0 < tls)
    {
        //printf("Late energy injection.\n");
        u = u1*pow(t0/t1, -0.5*(2+q)/(2-q));
        dr = (t1*dr1 + (2-q)*(t0/(u*u)-t1/(u1*u1))/16.0) / t0;

        *R0 = c*t0*(1-dr);
        *u0 = u;
        return;
    }

    //Energy injection is over!
    //Just use original solution with boosted energy.
    //printf("All energy injected.\n");
    Ei = E_inj(ts, L0, q, ts);
    C = sqrt(9/(16.0*M_PI) * (E0+Ei)/(rho0*c5));
    u = C * pow(t0, -1.5);

    *R0 = c*t0*(1 - 1/(16*u*u));
    *u0 = u;
}

void shockInitFind(double t0, double *R0, double *u0, double tRes, void *argv)
{
    struct shockParams *par = (struct shockParams *)argv;
    if(par->envType == ENV_ISM)
    {
        shockInitFindISM(t0, R0, u0, tRes, argv);
        return;
    }


    double E0 = par->E0;
    double Mej = par->Mej;
    double c2 = v_light * v_light;

    double t00, R00, u00;

    if(Mej > 0.0)
    {
        double gm1_coast = E0 / (Mej * c2);
        double g_coast = gm1_coast + 1;
        double u2_coast = gm1_coast * (gm1_coast + 2);
        double be2_coast = u2_coast / (g_coast * g_coast);

        double Msw = 3.0 * E0 / ((4*u2_coast + 3) * be2_coast * c2);

        //Rough deceleration radius
        double Rd = envRadiusPar(Msw, par);

        //shock 4-velocity
        double u_coast = sqrt(u2_coast);
        double uS = shockVel(u_coast);
        double beS = uS / sqrt(uS*uS + 1);
        double td = Rd / (beS * v_light);

        if(t0 < 0.001*td)
        {
            *R0 = beS * v_light * t0;
            *u0 = u_coast;
            return;
        }

        t00 = 0.001*td;
        R00 = beS * v_light * t00;
        u00 = u_coast;
    }
    else
    {
        //Setting ref ultra-relativistic velocity
        double u_UR = 1000;
        double g2_UR = u_UR*u_UR + 1;
        double be2_UR = u_UR*u_UR / g2_UR;

        //swept mass at u_UR;
        double M_UR = 3.0 * E0 / ((4*u_UR*u_UR + 3) * be2_UR * c2);

        //Rough UR radius
        double R_UR = envRadiusPar(M_UR, par);
        double rho_UR = envDensityPar(R_UR, par);

        // = 3-k for a power law rho ~ R^-k
        double volIdx = 4*M_PI * R_UR*R_UR*R_UR * rho_UR / M_UR;

        double t_UR = R_UR/v_light * (1.0 + 1.0/(4*(1+volIdx) * g2_UR));

        if(t0 < t_UR)
        {
            double R_apprx = v_light * t0;
            double M0 = envMassPar(R_apprx, par);
            double g2 = 0.75 * E0 / (M0 * c2);
            *R0 = R_apprx * (1 - 1.0/(4 * (1+volIdx) * g2));
            *u0 = sqrt(g2 - 1.0);
            return;
        }

        t00 = t_UR;
        R00 = R_UR;
        u00 = u_UR;
    }

    double dt, x[2], x0[2], k1[2], k2[2], k3[2], k4[2];

    x[0] = R00;
    x[1] = u00;

    double t = t00;
    double fac = pow(10, 1.0/tRes);
    int j;

    while(t < t0)
    {
        dt = (fac-1)*t;
        if(fac*t >= t0)
            dt = t0-t;
        x0[0] = x[0];
        x0[1] = x[1];
        Rudot2D(t, x, argv, k1);
        
        for(j=0; j<2; j++)
            x[j] = x0[j] + 0.5*dt*k1[j];
        Rudot2D(t, x, argv, k2);
        
        for(j=0; j<2; j++)
            x[j] = x0[j] + 0.5*dt*k2[j];
        Rudot2D(t, x, argv, k3);
        
        for(j=0; j<2; j++)
            x[j] = x0[j] + dt*k3[j];
        Rudot2D(t, x, argv, k4);
        
        for(j=0; j<2; j++)
            x[j] = x0[j] + dt*(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6.0;
        t *= fac;
    }

    *R0 = x[0];
    *u0 = x[1];
}

void shockInitFindISM(double t0, double *R0, double *u0, double tRes,
                      void *argv)
{
    struct shockParams *par = (struct shockParams *)argv;

    double E0 = par->E0;
    double Mej = par->Mej;
    double rho0 = par->rho0_env;
    double L0 = par->L0_inj;
    double q = par->q_inj;
    double ts = par->ts_inj;

    //printf("shockInitFind: t0=%.6lg E0=%.6lg Mej=%.6lg rho0=%.6lg\n",
    //        t0, E0, Mej, rho0);
    //printf("               L0=%.6lg q=%.6lg ts=%.6lg\n", L0, q, ts);

    double u;
    double c = v_light;
    double c5 = c*c*c*c*c;
    
    double C = sqrt(9/(16.0*M_PI) * E0/(rho0*c5));
    double tNR = pow(C, 2.0/3.0);

    //Time for deceleration
    double td;
    if(Mej > 0.0)
    {
        // gc, uc, beSc = coasting gamma, u, beta_shock
        double gcmo = E0/(Mej *c*c);
        double uc = sqrt(gcmo * (gcmo+2));
        double gc = sqrt(1+uc*uc);
        double beSc = 4*uc*gc / (4*uc*uc+3);
        double Rd = pow(9*gc*gc*Mej / (4*M_PI*(gc+1)*(4*uc*uc+3)*rho0), 1./3.);
        td = Rd / (beSc * c);
    }
    else
    {
        td = 0.0;
    }

    //Time for transition
    double ti;
    if(L0 > 0 && ts > 0)
    {
        double tei = te_inj(E0, L0, q, ts);
        ti = pow(16*C*C*tei, 0.25);
        if(tei < 0.0)
            ti = 1.0e20 * tNR; //Arbitrary val
    }
    else
        ti = 1.0e20 * tNR; //Arbitrary val

    //printf("               td=%.6lg ti=%.6lg tNR=%.6lg\n", td, ti, tNR);

    if(t0 < 0.01*tNR && t0 < 0.01*ti && t0 > 100*td)
    {
        //printf("the EASY way\n");
        u = C*pow(t0, -1.5);
        *R0 = c*t0*(1-1/(16*u*u));
        *u0 = u;
        return;
    }
    else if(t0 < 0.01 * td)
    {
        double gcmo = E0/(Mej *c*c);
        double uc = sqrt(gcmo * (gcmo+2));
        double gc = sqrt(1+uc*uc);
        double beSc = 4*uc*gc / (4*uc*uc+3);

        *R0 = beSc * c * t0;
        *u0 = uc;
        return;
    }
    
    double dt, x[2], x0[2], k1[2], k2[2], k3[2], k4[2];

    double t00, u00, R00;

    if(td > 0.0)
    {
        double gcmo = E0/(Mej *c*c);
        double uc = sqrt(gcmo * (gcmo+2));
        double gc = sqrt(1+uc*uc);
        double beSc = 4*uc*gc / (4*uc*uc+3);

        t00 = td<t0 ? 0.01*td : 0.01*t0;
        R00 = beSc * c * t0;
        u00 = uc;
    }
    else
    {
        t00 = ti < tNR ? 0.01*ti : 0.01*tNR;
        u00 = C*pow(t00, -1.5);
        R00 = c*t00*(1-1/(16*u00*u00));
    }

    //printf("t00=%.6le R00=%.6le u00=%.6le (tNR=%.3le ti=%.3le)\n", t00, R00, u00, tNR, ti);

    x[0] = R00;
    x[1] = u00;

    double t = t00;
    double fac = pow(10, 1.0/tRes);
    int j;

    while(t < t0)
    {
        dt = (fac-1)*t;
        if(fac*t >= t0)
            dt = t0-t;
        x0[0] = x[0];
        x0[1] = x[1];
        Rudot2D(t, x, argv, k1);
        
        for(j=0; j<2; j++)
            x[j] = x0[j] + 0.5*dt*k1[j];
        Rudot2D(t, x, argv, k2);
        
        for(j=0; j<2; j++)
            x[j] = x0[j] + 0.5*dt*k2[j];
        Rudot2D(t, x, argv, k3);
        
        for(j=0; j<2; j++)
            x[j] = x0[j] + dt*k3[j];
        Rudot2D(t, x, argv, k4);
        
        for(j=0; j<2; j++)
            x[j] = x0[j] + dt*(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6.0;
        t *= fac;
    }

    *R0 = x[0];
    *u0 = x[1];
}

void Rudot2D(double t, double *x, void *argv, double *xdot)
{
    struct shockParams *par = (struct shockParams *)argv;

    double Mej = par->Mej;
    double Einj = par->E0_refresh;
    double k = par->k_refresh;
    double umin = par->umin_refresh;
    double L0 = par->L0_inj;
    double q = par->q_inj;
    double ts = par->ts_inj;

    double R = x[0];
    double u = x[1];

    double g = sqrt(1+u*u);
    double be = u/g;

    double bes = 4*u*g/(4*u*u+3);

    double dRdt = bes * v_light;

    double dEdu = 0.0;
    if(Einj > 0.0 && u>umin)
        dEdu = -k*Einj*pow(u,-k-1);

    double dEdt = 0.0;
    double te = t - R/v_light;
    if(L0 > 0.0 && te < ts)
    {
        double gs2 = (4*u*u+3)*(4*u*u+3) / (8*u*u+9);
        dEdt = L_inj(te, L0, q, ts) / (gs2*(1+bes)); // Doppler factor (1-bes)
    }

    double rho = envDensityPar(R, par);
    double M = envMassPar(R, par);

    double dR_Esh = 4*M_PI/3.0 * (4*u*u+3)*be*be * rho * R*R * v_light*v_light;
    double du_Esh = 2.0 * (2*u*u+1)*(2*u*u+3)*u * M * v_light*v_light
                    / (3.0 * g*g*g*g);
    double du_Eej = be * Mej * v_light*v_light;
    double dudt = (-dR_Esh * dRdt + dEdt) / (du_Eej + du_Esh - dEdu);

    //double num = -16*M_PI/3.0 * rho0*R*R * be*u*u * v_light
    //                + dEdt/(v_light*v_light);
    //double denom = be*Mej
    //                + 8*M_PI*rho0*R*R*R*u*(2*u*u+1)*(2*u*u+3)/(9*g*g*g*g)
    //                - dEdu/(v_light*v_light);
    //double dudt = num/denom;

    xdot[0] = dRdt;
    xdot[1] = dudt;
}

void RuThdot3D(double t, double *x, void *argv, double *xdot)
{
    struct shockParams *par = (struct shockParams *)argv;

    double Mej = par->Mej;
    //double rho0 = par->rho0_env;
    double Einj = par->E0_refresh;
    double k = par->k_refresh;
    double umin = par->umin_refresh;
    double L0 = par->L0_inj;
    double q = par->q_inj;
    double ts = par->ts_inj;
    double thC = par->thetaCore;
    double th0 = par->theta0;
    double thCg = par->thetaCoreGlobal;
    int spread = par->spread;

    double R = x[0];
    double u = x[1];
    double th = x[2];

    double g = sqrt(1+u*u);
    double be = u/g;
    double bes = 4*u*g/(4*u*u+3);

    double sinth = sin(0.5*th);
    double costh = cos(0.5*th);
    double om = 2*sinth*sinth;

    double dRdt = 4*u*g/(4*u*u+3) * v_light;

    double rho = envDensityPar(R, par);
    double M = envMassPar(R, par);
    double env3mk = 4*M_PI*rho*R*R*R / M;

    double dThdt = 0.0;
    //TODO: verify spreading procedure 190401
    //if(spread && th < 0.5*M_PI && u < 1)
    //if(spread && th < 0.5*M_PI && u*3.0*th < 1)
    if(spread)
    {
        if(spread == 1)
        {
            double Q0 = 2.0;
            double Q = sqrt(2.0)*3.0;

            if(th < 0.5*M_PI && Q0*u*thC < 1)
            {
                double e = u*u/(g+1); // specific internal energy == gamma-1

                //Sound speed from trans-relativistic EoS.
                double cs = v_light * sqrt(e*(2+e)*(5+8*e+4*e*e)
                                            / (3*(1+e)*(1+e)*(1+2*e)*(3+2*e)));
                double fac = u*thC*Q < 1.0 ? 1.0 : Q*(1-Q0*u*thC) / (Q-Q0);
                /*
                if(fac < 1.0)
                {
                    double sharp = 3.0;
                    double f0 = exp(-sharp);
                    fac = (exp(sharp*(fac-1.0))-f0) / (1.0-f0);
                }
                */
                //fac = 0.0;
                dThdt = fac * cs / (R*g);
            }
        }
        else if(spread == 2)
        {
            double Q0 = 2.0;
            double Q = sqrt(2)*3.0;

            if(th < 0.5*M_PI && Q0*u*th0 < 1)
            {
                double e = u*u/(g+1); // specific internal energy == gamma-1

                //Sound speed from trans-relativistic EoS.
                double cs = v_light * sqrt(e*(2+e)*(5+8*e+4*e*e)
                                            / (3*(1+e)*(1+e)*(1+2*e)*(3+2*e)));
                double fac = u*th0*Q < 1.0 ? 1.0 : Q*(1-Q0*u*th0) / (Q-Q0);
                /*
                if(fac < 1.0)
                {
                    double sharp = 3.0;
                    double f0 = exp(-sharp);
                    fac = (exp(sharp*(fac-1.0))-f0) / (1.0-f0);
                }
                */
                //fac = 0.0;
                dThdt = fac * cs / (R*g);
            }
        }
        else if(spread == 3)
        {
            double Q0 = 2.0;
            double Q = sqrt(2)*3.0;

            if(th < 0.5*M_PI && Q0*u*thCg < 1)
            {
                double e = u*u/(g+1); // specific internal energy == gamma-1

                //Sound speed from trans-relativistic EoS.
                double cs = v_light * sqrt(e*(2+e)*(5+8*e+4*e*e)
                                            / (3*(1+e)*(1+e)*(1+2*e)*(3+2*e)));
                double fac = u*thCg*Q < 1.0 ? 1.0 : Q*(1-Q0*u*thCg) / (Q-Q0);
                /*
                if(fac < 1.0)
                {
                    double sharp = 3.0;
                    double f0 = exp(-sharp);
                    fac = (exp(sharp*(fac-1.0))-f0) / (1.0-f0);
                }
                */
                //fac = 0.0;
                dThdt = fac * cs / (R*g);
            }
        }
        else if(spread == 4)
        {
            double bew = 0.5*sqrt((2*u*u+3)/(4*u*u+3))*bes/g;
            dThdt = bew * v_light / R;
        }
        else if(spread == 5)
        {
            double Q0 = 2.0;
            double Q = sqrt(2.0)*3.0;

            if(th < 0.5*M_PI && Q0*u*thC < 1)
            {
                double bew = 0.5*sqrt((2*u*u+3)/(4*u*u+3))*bes/g;
                double fac = u*thC*Q < 1.0 ? 1.0 : Q*(1-Q0*u*thC) / (Q-Q0);
                dThdt = fac * bew * v_light / R;
            }
        }
        else if(spread == 6)
        {
            double Q0 = 2.0;
            double Q = sqrt(2.0)*3.0;

            if(th < 0.5*M_PI && Q0*u*thCg < 1)
            {
                double bew = 0.5*sqrt((2*u*u+3)/(4*u*u+3))*bes/g;
                double fac = u*thCg*Q < 1.0 ? 1.0 : Q*(1-Q0*u*thCg) / (Q-Q0);
                dThdt = fac * bew * v_light / R;
            }
        }
        else if(spread == 7)
        {
            double envk = 3 - env3mk;
            double Q = sqrt(2.0)*env3mk;
            double Q0;
            if(envk < 2.0)
                Q0 = 0.5 * (2.0 * (2.0 - envk) + 1.4 * envk);
            else
                Q0 = 1.4 * (3.0 - envk);

            if(th < 0.5*M_PI && Q0*u*thC < 1)
            {
                double bew = 0.5*sqrt((2*u*u+3)/(4*u*u+3))*bes/g;
                double fac = u*thC*Q < 1.0 ? 1.0 : Q*(1-Q0*u*thC) / (Q-Q0);
                if (th0 < thC)
                    fac *= tan(0.5*th0)/tan(0.5*thC); //th0/thC;
                dThdt = fac * bew * v_light / R;
            }
        }
        else if(spread == 8)
        {
            double envk = 3 - env3mk;
            double Q = sqrt(2.0)*env3mk;
            double Q0;
            if(envk < 2.0)
                Q0 = 0.5 * (2.0 * (2.0 - envk) + 1.4 * envk);
            else
                Q0 = 1.4 * (3.0 - envk);

            if(th < 0.5*M_PI && Q0*u*thCg < 1)
            {
                double bew = 0.5*sqrt((2*u*u+3)/(4*u*u+3))*bes/g;
                double fac = u*thCg*Q < 1.0 ? 1.0 : Q*(1-Q0*u*thCg) / (Q-Q0);
                if (th0 < thCg)
                    fac *= tan(0.5*th0)/tan(0.5*thCg); //th0/thC;
                dThdt = fac * bew * v_light / R;
            }
        }
        else if(spread == 9)
        {
            //Fernandez & Lamb 2021 spreading.
            if(th < 0.5*M_PI && u*thC < 100)
            {
                double gs = 1 / sqrt((1-bes)*(1+bes));
                double fac = 1.0;
                if (th0 < thC)
                    fac *= tan(0.5*th0)/tan(0.5*thC); //th0/thC;
                dThdt = fac / (gs*gs * th) * dRdt / R;
            }
        }
        else if(spread == 11)
        {
            double envk = 3 - env3mk;
            double Q = sqrt(2.0)*env3mk;
            double Q0;
            if(envk < 2.0)
                Q0 = 0.5 * (2.0 * (2.0 - envk) + 1.4 * envk);
            else
                Q0 = 1.4 * (3.0 - envk);

            if(th < 0.5*M_PI && Q0*u*thC < 1)
            {
                double bew = 0.5*sqrt((2*u*u+3)/(4*u*u+3))*bes/g;
                double fac = u*thC*Q < 1.0 ? 1.0 : Q*(1-Q0*u*thC) / (Q-Q0);
                if (th0 < thC)
                    fac *= tan(0.5*th0)/tan(0.5*thC); //th0/thC;
                double gs = 1.0 / sqrt((1-bes)*(1+bes));
                double dThdt_lum = fac * v_light / (R * g);
                dThdt = 4.0 * fac * bew * v_light / R;
                if(dThdt > dThdt_lum)
                    dThdt = dThdt_lum;
            }
        }
    }

    double dEdu = 0.0;
    if(Einj > 0.0 && u>umin)
        dEdu = -k*Einj*pow(u,-k-1);

    double dEdt = 0.0;
    double te = t - R/v_light;
    if(L0 > 0.0 && te < ts)
    {
        double gs2 = (4*u*u+3)*(4*u*u+3) / (8*u*u+9);
        dEdt = L_inj(te, L0, q, ts) / (gs2*(1+bes)); // Doppler factor (1-bes)
    }

    double c2 = v_light*v_light;

    double dR_Esh = 4*M_PI/3.0 * (4*u*u+3)*be*be * rho * R*R * om * c2;
    double du_Esh = 2.0*(2*u*u+1)*(2*u*u+3)*u * M * om * c2 / (3.0 * g*g*g*g);
    double dTh_Esh = (4*u*u+3)*be*be * M * 2*sinth*costh * c2 / 3.0;
    double du_Eej = be * Mej * c2;
    double dudt = (-dR_Esh * dRdt - dTh_Esh * dThdt + dEdt)
                   / (du_Eej + du_Esh - dEdu);

    /*
    double num = -16*M_PI/3.0 * om*rho0*R*R * be*u*u * v_light
                 -8*M_PI/9.0*rho0*R*R*R*(4*u*u+3)*be*be*sinth*costh*dThdt
                    + om*dEdt/(v_light*v_light);
    double denom = be*Mej
                    + 8*M_PI*om*rho0*R*R*R*u*(2*u*u+1)*(2*u*u+3)/(9*g*g*g*g)
                    - dEdu/(v_light*v_light);
    double dudt = num/denom;
    */

    xdot[0] = dRdt;
    xdot[1] = dudt;
    xdot[2] = dThdt;
}

void shockEvolveRK4(double *t, double *R, double *u, int N, double R0, 
                    double u0, void *args)
{
    int i,j;

    double dt, x[2], x0[2], k1[2], k2[2], k3[2], k4[2];

    R[0] = R0;
    u[0] = u0;

    for(i=0; i<N-1; i++)
    {
        dt = t[i+1] - t[i];
        x0[0] = R[i];
        x0[1] = u[i];
        Rudot2D(t[i], x0, args, k1);
        
        for(j=0; j<2; j++)
            x[j] = x0[j] + 0.5*dt*k1[j];
        Rudot2D(t[i], x, args, k2);
        
        for(j=0; j<2; j++)
            x[j] = x0[j] + 0.5*dt*k2[j];
        Rudot2D(t[i], x, args, k3);
        
        for(j=0; j<2; j++)
            x[j] = x0[j] + dt*k3[j];
        Rudot2D(t[i], x, args, k4);
        
        for(j=0; j<2; j++)
            x[j] = x0[j] + dt*(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6.0;

        R[i+1] = x[0];
        u[i+1] = x[1];
    }
}
void shockEvolveSpreadRK4(double *t, double *R, double *u, double *th, int N,
                            double R0, double u0, double th0, void *args)
{
    int i,j;

    double dt, x[3], x0[3], k1[3], k2[3], k3[3], k4[3];

    R[0] = R0;
    u[0] = u0;
    th[0] = th0;

    for(i=0; i<N-1; i++)
    {
        dt = t[i+1] - t[i];
        x0[0] = R[i];
        x0[1] = u[i];
        x0[2] = th[i];
        RuThdot3D(t[i], x0, args, k1);
        
        for(j=0; j<3; j++)
            x[j] = x0[j] + 0.5*dt*k1[j];
        RuThdot3D(t[i], x, args, k2);
        
        for(j=0; j<3; j++)
            x[j] = x0[j] + 0.5*dt*k2[j];
        RuThdot3D(t[i], x, args, k3);
        
        for(j=0; j<3; j++)
            x[j] = x0[j] + dt*k3[j];
        RuThdot3D(t[i], x, args, k4);
        
        for(j=0; j<3; j++)
            x[j] = x0[j] + dt*(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6.0;

        R[i+1] = x[0];
        u[i+1] = x[1];
        if(x[2] > 0.5*M_PI)
            th[i+1] = 0.5*M_PI;
        else
            th[i+1] = x[2];
    }
}

void setup_shockParams(struct shockParams *pars, int spread,
                       double E0, double Mej,
                       int envType, double rho0_env, double R0_env,
                       double k_env, double rho1_env,
                       double L0_inj, double q_inj, double t0_inj,
                       double ts_inj,
                       double E0_refresh, double k_refresh, double umin_refresh,
                       double thetaCore, double theta0,
                       double thetaCoreGlobal)
{
    pars->spread = spread;
    pars->E0 = E0;
    pars->Mej = Mej;
    pars->envType = envType;
    pars->rho0_env = rho0_env;
    pars->R0_env = R0_env;
    pars->k_env = k_env;
    pars->rho1_env = rho1_env;
    pars->L0_inj = L0_inj;
    pars->q_inj = q_inj;
    pars->t0_inj = t0_inj;
    pars->ts_inj = ts_inj;
    pars->E0_refresh = E0_refresh;
    pars->k_refresh = k_refresh;
    pars->umin_refresh = umin_refresh;
    pars->thetaCore = thetaCore;
    pars->theta0 = theta0;
    pars->thetaCoreGlobal = thetaCoreGlobal;
}
