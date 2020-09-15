#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "integrate.h"

#define KMAX 20

double trap(double (*f)(double, void *), double xa, double xb, int N,
            void *args, int (*errf)(void *))
{
    double dx = (xb - xa)/N;
    double I = 0.5*(f(xa, args) + f(xb, args));
    int i;
    for(i=1; i<N; i++)
    {
        I += f(xa + i*dx, args);
        if(errf(args))
            return 0.0;
    }
    return I*dx;
}

double simp(double (*f)(double, void *), double xa, double xb, int N,
            void *args, int (*errf)(void *))
{
    if(N%2 == 1)
        N -= 1;

    double dx = (xb - xa)/N;
    double I1, I2, I3;
    int i;
    I1 = f(xa, args) + f(xb, args);
    if(errf(args))
        return 0.0;
    I2 = 0.0;
    for(i=1; i<N; i+=2)
    {
        I2 += f(xa + i*dx, args);
        if(errf(args))
            return 0.0;
    }
    I3 = 0.0;
    for(i=2; i<N; i+=2)
    {
        I3 += f(xa + i*dx, args);
        if(errf(args))
            return 0.0;
    }
    return (I1 + 4*I2 + 2*I3) * dx / 3.0;
}

double romb(double (*f)(double, void *), double xa, double xb, int N,
            double atol, double rtol, void *args, int *Neval, double *eps,
            int verbose, int (*errf)(void *))
{
    double R[KMAX];

    int m, k, k0, Nk;
    long fpm;
    double hk, Rp, err;

    double maxFracChange = 0.1;

    hk = xb - xa;
    Nk = 1;
    R[KMAX-1] = 0.5*(xb-xa)*(f(xa, args) + f(xb, args));
    if(errf(args))
        return 0.0;
    R[0] = R[KMAX-1];

    if(Neval != NULL)
        *Neval = 2;

    for(k=1; k<KMAX; k++)
    {
        k0 = KMAX-k-1;
        hk *= 0.5;
        Nk *= 2;

        Rp = 0.0;
        for(m=1; m<Nk; m+=2)
        {
            Rp += f(xa + m*hk, args);
            if(errf(args))
                return 0.0;
        }
        R[k0] = 0.5*R[k0+1] + hk*Rp;
        if(Neval != NULL)
            *Neval += Nk/2;

        fpm = 1;
        for(m=1; m<=k; m++)
        {
            fpm *= 4;
            R[k0+m] = (fpm*R[k0+m-1] - R[k0+m]) / (fpm - 1);
        }
        err = (R[KMAX-1] - R[0]) / (fpm - 1);
        double lastVal = R[0];
        R[0] = R[KMAX-1];

        double fracChange = fabs((R[0] - lastVal) / lastVal);

        if(eps != NULL)
            *eps = err;
        if(verbose)
            printf("level %d:  Neval=%d  I=%.6lg  fracDelta=%.3lg"
                   " err=%.6lg  tol=%.6lg\n",
                    k, Nk+1, R[0], fracChange, err, atol+rtol*fabs(R[0]));

        //printf("      k%d: I=%.6le err=%.6le frac=%.3le\n", k, R[0], err, 
        //        fabs(err) / (atol + rtol*fabs(R[0])));

        if((fabs(err) < atol + rtol*fabs(R[0]))
                && fracChange < maxFracChange)
            break;

        if(N > 1 && Nk >= N)
            break;
    }

    return R[0];
}

double trap_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, Mesh3 *mout, int verbose, int (*errf)(void *))
{
    double I = m3_adapt(f, xa, xb, Nmax, trapProcessInterval,
                        trapSplitInterval, atol, rtol, args,
                        Neval, eps, mout, verbose, errf);
    return I;
}

double simp_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                  double atol, double rtol, void *args, int *Neval,
                  double *eps, Mesh5 *mout, int verbose, int (*errf)(void *))
{
    double I = m5_adapt(f, xa, xb, Nmax, simpProcessInterval,
                        simpSplitInterval, atol, rtol, args,
                        Neval, eps, mout, verbose, errf);
    return I;
}

int trapProcessInterval(double (*f)(double, void *), void *args, Interval3 *i,
                        int (*errf)(void *))
{
    double fa = i->fa;
    double fb = i->fb;
    double fm = f(0.5*(i->a+i->b), args);
    if(errf(args))
        return 1;
    i->fm = fm;

    double h = 0.5*(i->b - i->a);

    double I0 = 0.5*(fa + fb) * 2*h;
    double I1 = 0.5*(fa + 2*fm + fb) * h;

    double err = (I1 - I0) / 3.0;
    i->err = fabs(err);
    i->I = I1 + err;

    return 1;
}

int trapSplitInterval(double (*f)(double, void *), void *args,
                      Interval3 *i0, Interval3 *i1, Interval3 *i2,
                      int (*errf)(void *))
{
    double x = 0.5*(i0->a + i0->b);
    i1->a = i0->a;
    i1->b = x;
    i2->a = x;
    i2->b = i0->b;

    i1->fa = i0->fa;
    i1->fb = i0->fm;
    i2->fa = i0->fm;
    i2->fb = i0->fb;

    int n = 0;
    n += trapProcessInterval(f, args, i1, errf);
    n += trapProcessInterval(f, args, i2, errf);

    return n;
}

int simpProcessInterval(double (*f)(double, void *), void *args, Interval5 *i,
                        int (*errf)(void *))
{
    double fa = i->fa;
    double fb = i->fb;
    double fm = i->fm;
    double fl = f(0.75*i->a+0.25*i->b, args);
    if(errf(args))
        return 0.0;
    double fr = f(0.25*i->a+0.75*i->b, args);
    if(errf(args))
        return 0.0;
    i->fl = fl;
    i->fr = fr;

    double h = 0.25*(i->b - i->a);

    double I0 = 2*h * (fa + 4*fm + fb)/3.0;
    double I1 = h * (fa + 4*fl + 2*fm + 4*fr + fb)/3.0;

    double err = (I1 - I0) / 15.0;
    i->err = fabs(err);
    i->I = I1 + err;

    return 2;
}

int simpSplitInterval(double (*f)(double, void *), void *args,
                      Interval5 *i0, Interval5 *i1, Interval5 *i2,
                      int (*errf)(void *))
{
    double xm = 0.5*(i0->a + i0->b);
    i1->a = i0->a;
    i1->b = xm;
    i2->a = xm;
    i2->b = i0->b;

    i1->fa = i0->fa;
    i1->fm = i0->fl;
    i1->fb = i0->fm;
    i2->fa = i0->fm;
    i2->fm = i0->fr;
    i2->fb = i0->fb;

    int n = 0;
    n += simpProcessInterval(f, args, i1, errf);
    n += simpProcessInterval(f, args, i2, errf);

    return n;
}

double m3_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                 int (*processInterval)(double (*f)(double, void*), void *,
                                         Interval3 *, int (*errf)(void *)),
                 int (*splitInterval)(double (*f)(double, void *), void *,
                      Interval3 *, Interval3 *, Interval3 *,
                      int (*errf)(void *)),
                 double atol, double rtol, void *args, int *Neval,
                 double *eps, Mesh3 *mout, int verbose, int (*errf)(void *))
{
    Mesh3 m;
    mesh3Init(&m);

    Interval3 i = {.a=xa, .b=xb, .I=0, .err=0, .fa=0, .fm=0, .fb=0};
    i.fa = f(xa, args);
    if(errf(args))
    {
        mesh3Free(&m);
        return 0.0;
    }
    i.fb = f(xb, args);
    if(errf(args))
    {
        mesh3Free(&m);
        return 0.0;
    }
    int n = 2;
    n += processInterval(f, args, &i, errf);
    if(errf(args))
    {
        mesh3Free(&m);
        return 0.0;
    }

    mesh3Insert(&m, i);

    double I = i.I;
    double err = i.err;
    int num_intervals = 1;
    int last_check = 1;

    while(n < Nmax && err >= atol + rtol*fabs(I))
    {
        i = mesh3Extract(&m);

        Interval3 i1;
        Interval3 i2;
        n += splitInterval(f, args, &i, &i1, &i2, errf);
        if(errf(args))
        {
            mesh3Free(&m);
            return 0.0;
        }
        mesh3Insert(&m, i1);
        mesh3Insert(&m, i2);
        num_intervals++;

        if(num_intervals == 2*last_check)
        {
            err = mesh3TotalError(&m);
            I = mesh3TotalIntegral(&m);
            last_check = num_intervals;
        }
        else
        {
            err += i1.err + i2.err - i.err;
            I += i1.I + i2.I - i.I;
        }
        
        if(verbose)
            printf("Num Intervals: %d - I=%.12lg  err=%.3lg  tol=%.3lg"
                   "  meshOk=%d\n",
                   num_intervals, I, err, atol + rtol*fabs(I), mesh3Check(&m));
    }

    I = mesh3TotalIntegral(&m);

    if(Neval != NULL)
        *Neval = n;

    if(eps != NULL)
    {
        err = mesh3TotalError(&m);
        *eps = err;
    }

    if(mout == NULL)
        mesh3Free(&m);
    else
        *mout = m;

    return I;
}

double m5_adapt(double (*f)(double, void *), double xa, double xb, int Nmax,
                 int (*processInterval)(double (*f)(double, void*), void *,
                                         Interval5 *, int (*errf)(void *)),
                 int (*splitInterval)(double (*f)(double, void *), void *,
                      Interval5 *, Interval5 *, Interval5 *,
                      int (*errf)(void *)),
                 double atol, double rtol, void *args, int *Neval,
                 double *eps, Mesh5 *mout, int verbose, int (*errf)(void *))
{
    Mesh5 m;
    mesh5Init(&m);

    Interval5 i = {.a=xa, .b=xb, .I=0, .err=0,
                   .fa=0, .fl=0, .fm=0, .fr=0, .fb=0};
    i.fa = f(xa, args);
    if(errf(args))
    {
        mesh5Free(&m);
        return 0.0;
    }
    i.fm = f(0.5*(xa+xb), args);
    if(errf(args))
    {
        mesh5Free(&m);
        return 0.0;
    }
    i.fb = f(xb, args);
    if(errf(args))
    {
        mesh5Free(&m);
        return 0.0;
    }
    int n = 3;
    n += processInterval(f, args, &i, errf);
    if(errf(args))
    {
        mesh5Free(&m);
        return 0.0;
    }

    mesh5Insert(&m, i);

    double I = i.I;
    double err = i.err;
    int num_intervals = 1;
    int last_check = 1;

    while(n < Nmax && err > atol + rtol*fabs(I))
    {
        i = mesh5Extract(&m);

        Interval5 i1;
        Interval5 i2;
        n += splitInterval(f, args, &i, &i1, &i2, errf);
        if(errf(args))
        {
            mesh5Free(&m);
            return 0.0;
        }
        mesh5Insert(&m, i1);
        mesh5Insert(&m, i2);
        num_intervals++;

        if(num_intervals == 2*last_check)
        {
            err = mesh5TotalError(&m);
            I = mesh5TotalIntegral(&m);
            last_check = num_intervals;
        }
        else
        {
            err += i1.err + i2.err - i.err;
            I += i1.I + i2.I - i.I;
        }
        
        if(verbose)
            printf("Num Intervals: %d - I=%.12lg  err=%.3lg  tol=%.3lg"
                   "  meshOk=%d\n",
                   num_intervals, I, err, atol + rtol*fabs(I), mesh5Check(&m));
    }

    I = mesh5TotalIntegral(&m);

    if(Neval != NULL)
        *Neval = n;

    if(eps != NULL)
    {
        err = mesh5TotalError(&m);
        *eps = err;
    }

    if(mout == NULL)
        mesh5Free(&m);
    else
        *mout = m;

    return I;
}
