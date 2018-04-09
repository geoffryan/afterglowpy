#include <stdio.h>
#include <math.h>
#include "offaxis_struct.h"
#include "shockEvolution.h"

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

    double R = x[0];
    double u = x[1];

    double g = sqrt(1+u*u);
    double be = u/g;

    double dRdt = 4*u*g/(4*u*u+3) * v_light;

    double dEdu = 0.0;
    if(Einj > 0.0 && u>umin)
        dEdu = -k*Einj*pow(u,-k-1);

    double dEdt = 0.0;
    if(L0 > 0.0)
        dEdt = L0*pow(t/1.0e3, -q);

    double num = -16*M_PI/3.0 * rho0*R*R * be*u*u * v_light
                    + dEdt/(v_light*v_light);
    double denom = be*Mej + 8*M_PI*rho0*R*R*R*u*(2*u*u+1)*(2*u*u+3)/(9*g*g*g*g)
                    - dEdu/(v_light*v_light);
    double dudt = num/denom;

    xdot[0] = dRdt;
    xdot[1] = dudt;
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
