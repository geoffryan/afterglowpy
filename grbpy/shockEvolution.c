#include <stdio.h>
#include <math.h>
#include "offaxis_struct.h"
#include "shockEvolution.h"

double shockVel(double u)
{
    return 4*u*sqrt((u*u+1) / (8*u*u+9));
}

void Rudot2D(double t, double *x, void *argv, double *xdot)
{
    double *args = (double *)argv;

    double Mej = args[1];
    double rho0 = args[2];
    double Einj = args[3];
    double k = args[4];
    double umin = args[5];
    double L0 = args[6];
    double q = args[7];
    double ts = args[8];

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
        dEdt = L0*pow(te/1.0e3, -q) / (gs2*(1+bes)); // Doppler factor (1-bes)
    }

    double num = -16*M_PI/3.0 * rho0*R*R * be*u*u * v_light
                    + dEdt/(v_light*v_light);
    double denom = be*Mej + 8*M_PI*rho0*R*R*R*u*(2*u*u+1)*(2*u*u+3)/(9*g*g*g*g)
                    - dEdu/(v_light*v_light);
    double dudt = num/denom;

    xdot[0] = dRdt;
    xdot[1] = dudt;
}

void RuThdot3D(double t, double *x, void *argv, double *xdot)
{
    double *args = (double *)argv;

    double Mej = args[1];
    double rho0 = args[2];
    double Einj = args[3];
    double k = args[4];
    double umin = args[5];
    double L0 = args[6];
    double q = args[7];
    double ts = args[8];

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

    double dThdt = 0.0;
    if(th < 0.5*M_PI && u < 1)
    {
        double e = u*u/(g+1); // specific internal energy == gamma-1

        //Sound speed from trans-relativistic EoS.
        double cs = v_light * sqrt(e*(2+e)*(5+8*e+4*e*e)
                                    / (3*(1+e)*(1+e)*(1+2*e)*(3+2*e)));
        dThdt = cs / (R*g);
    }

    double dEdu = 0.0;
    if(Einj > 0.0 && u>umin)
        dEdu = -k*Einj*pow(u,-k-1);

    double dEdt = 0.0;
    double te = t - R/v_light;
    if(L0 > 0.0 && te < ts)
    {
        double ti0 = 1.0e3;
        double gs2 = (4*u*u+3)*(4*u*u+3) / (8*u*u+9);
        dEdt = L0*pow(te/ti0, -q) / (gs2*(1+bes)); // Doppler factor (1-bes)
    }

    double num = -16*M_PI/3.0 * om*rho0*R*R * be*u*u * v_light
                 -8*M_PI/9.0*rho0*R*R*R*(4*u*u+3)*be*be*sinth*costh*dThdt
                    + om*dEdt/(v_light*v_light);
    double denom = be*Mej
                    + 8*M_PI*om*rho0*R*R*R*u*(2*u*u+1)*(2*u*u+3)/(9*g*g*g*g)
                    - dEdu/(v_light*v_light);
    double dudt = num/denom;

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
        th[i+1] = x[2];
    }
}
