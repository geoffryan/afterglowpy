import math
import numpy as np
import scipy.integrate as integrate
from . import shock

c = 2.99792458e10
me = 9.1093897e-28
mp = 1.6726231e-24
h = 6.6260755e-27
hbar = 1.05457266e-27
ee = 4.803e-10
sigmaT = 6.65e-25

Msun = 1.98892e33
cgs2mJy = 1.0e26
mJy2cgs = 1.0e-26
deg2rad = np.pi/180.0
rad2deg = 180.0/np.pi
day2sec = 86400.0
sec2day = 1.0/day2sec
parsec = 3.0857e18
Hz2eV = 4.13566553853599e-15
eV2Hz = 1.0/Hz2eV

intrtol = 1.0e-3


def dfdt(f, t, umax, umin, Ei, q, Mej_solar, Vej0, rho0):

    R = f[0]
    u = f[1]

    g = math.sqrt(u*u+1)
    beta = u/g

    Mej = Mej_solar * Msun

    dRdt = 4*c* u * g / (4*u*u + 3)

    dEdt = 0.0
    #if ts > 0.0 and t < ts:
    #    dEdt = Ei/ts * math.pow(t/ts, q)

    dEdu = 0.0
    if u > umin and u < umax:
        dEdu = -q*Ei/(c*c) * math.pow(u, -(q+1))


    A_g1 = 4.0/3.0 * np.pi * rho0 * R*R/(g*g) * (4*g*g*g*g - 5*g*g + 1) * dRdt
    A_g2 = beta * (Mej - 4.0/3.0 * rho0 * Vej0 * (19*g*g + 4*g)
            + 8.0/9.0*np.pi*rho0 * R*R*R * (4*g*g*g*g - 1.0)
                / (g*g*g))

    #ai = 5.0/3.0
    #A_g1 = np.pi * R*R * rho0 * u*u*(1+(g-1)*ai)*(1+g*ai) / (g*g*(ai-1)) * dRdt
    #A_g2 = beta * (Mej + np.pi * R*R*R * rho0 * (2 + 2*(g*g*g+g-1)*ai
    #                    + g*(g-1)*(1+g+2*g*g)*ai*ai) / (3*g*g*g*(ai-1)))

    dudt = (dEdt/(c*c) - A_g1) / (A_g2 - dEdu)

    return np.array([dRdt,dudt])

def dP(theta, amu, ate, au, ar, nu, n0, p, epsE, epsB, ksiN, specType):

    mu = math.cos(theta)
    ib = np.searchsorted(amu, mu)
    N = amu.shape[0]
    if ib <= 0:
        ib = 1
    elif ib >= N:
        ib = N-1
    ia = ib-1

    te = ((mu-amu[ia])*ate[ib] + (amu[ib]-mu)*ate[ia]) / (amu[ib]-amu[ia])
    #u = ((mu-amu[ia])*au[ib] + (amu[ib]-mu)*au[ia]) / (amu[ib]-amu[ia])
    #r = ((mu-amu[ia])*ar[ib] + (amu[ib]-mu)*ar[ia]) / (amu[ib]-amu[ia])
    #te = ate[ia]*math.pow(mu/amu[ia], math.log(ate[ib]/ate[ia])
    #                                    / math.log(amu[ib]/amu[ia]))
    u = au[ia]*math.pow(te/ate[ia], math.log(au[ib]/au[ia]) 
                                        / math.log(ate[ib]/ate[ia]))
    r = ar[ia]*math.pow(te/ate[ia], math.log(ar[ib]/ar[ia])
                                        / math.log(ate[ib]/ate[ia]))

    g = math.sqrt(u*u+1)
    beta = u/g

    us = 4*u*g / math.sqrt(8*u*u+9)
    gs = math.sqrt(us*us + 1)
    betashock = us / gs

    nprime = 4.0 * n0 * g
    eth = u*u/(g+1.0) * nprime * mp*c*c # (g-1)*nprime * mp*c*c
    B = math.sqrt(8*np.pi*eth*epsB)
    a = 1.0 - mu*beta  # beaming factor
    ashock = 1.0 - mu*betashock # shock beaming factor

    DR = math.fabs(r / (12.0*g*g*ashock))

    nuprime = g*a * nu
    gm = (2.0-p)/(1.0-p) * epsE * eth / (ksiN * nprime * me*c*c)
    gcs = 6*np.pi * g * me*c / (sigmaT * B*B * te)

    gr = gcs/gm
    y = beta*epsE/epsB

    if specType == 1:
        if gr <= 1.0 or gr*gr-gr-y <= 0.0:
            #Fast cooling (with I.C.)
            X = 0.5 * (1 + math.sqrt(1+4*y))
        else:
            #Slow cooling (with I.C.)
            b = y*math.pow(gr, 2-p)
            Xa = 1+b
            Xb = math.pow(b, 1.0/(4.0-p)) + 1.0/(4.0-p)
            s = b*b/(b*b+1.0)
            X0 = Xa * math.pow(Xb/Xa, s)
            X = X0
            imax = 5
            for i in range(imax):
                po = math.pow(X, p-2)
                f = X*X - X - b*po
                df = 2*X - 1 - p*(p-2)*po/X
                dX = -f/df
                X += dX
                if math.fabs(dX) < 1e-4 * X:
                    break
        gc = gcs / X 
    else:
        gc = gcs

    num = 3*gm*gm * ee * B / (4*np.pi * me*c)
    nuc = 3*gc*gc * ee * B / (4*np.pi * me*c)
    em = ksiN * nprime * B

    if num < nuc:
        if nuprime < num:
            freq = math.pow(nuprime/num, 1.0/3.0)
        elif nuprime < nuc:
            freq = math.pow(nuprime/num, 0.5*(1.0-p))
        else:
            freq = math.pow(nuc/num, 0.5*(1.0-p)) * math.pow(nuprime/nuc,
                                                                    -0.5*p)
    else:
        if nuprime < nuc:
            freq = math.pow(nuprime/nuc, 1.0/3.0)
        elif nuprime < num:
            freq = math.sqrt(nuc/nuprime)
        else:
            freq = math.sqrt(nuc/num) * math.pow(nuprime/num, -0.5*p)

    return 2*np.pi * r*r*math.sin(theta) * DR * em * freq / (g*g*a*a)

def rk4(dfdt, x0, at, args):

    #print("In RK4")

    n = x0.shape[0]
    N = at.shape[0]
    ax = np.empty((N,n))

    ax[0,:] = x0

    for i in range(N-1):
        x = ax[i,:]
        dt = at[i+1] - at[i]
        k1 = dfdt(x, at[i], *args)
        k2 = dfdt(x + 0.5*k1*dt, at[i], *args)
        k3 = dfdt(x + 0.5*k2*dt, at[i], *args)
        k4 = dfdt(x + k3*dt, at[i], *args)
        ax[i+1,:] = x + dt*(k1+2*k2+2*k3+k4)/6.0

    return ax



def fluxDensity(t, nu, jetType, specType, umax, umin, Ei, k, Mej_solar,
                        n0, p, epsE, epsB, ksiN, dL):

    r0 = 1.0e9
    rho0 = mp * n0
    Mej = Mej_solar * Msun
    u0 = umax
    g0 = math.sqrt(1+u0*u0)
    bes0 = 4*u0*g0 / (4*u0*u0+3)
    Rd = math.pow(9*g0*g0*Mej / (4*np.pi*(g0+1)*(4*u0*u0+3)*rho0), 1./3.)
    td = Rd / (bes0 * c)

    t0 = min(1.0e-2*td, 5.0e-1 * g0*g0*t.min(), 5.0e-1 * t.min()/(1+bes0))
    t1 = 2. * g0*g0*t.max()

    NT = int(1000 * math.log10(t1/t0))
    #print("{0:.3e} {1:.3e} {2:.3e} {3:d}".format(t0, t1, t1/t0, NT))

    r0 = bes0*c*t0

    Vej0 = 4.0/3.0*np.pi*r0*r0*r0

    #r0 = 1.0e9
    #t0 = r0 / c
    #t1 = 1.0e12 * t0
    #NT = 15000
    ate = np.logspace(math.log10(t0), math.log10(t1), num=NT, base=10.0)

    """
    f0 = np.array([r0, umax])
    args = (umax, umin, Ei, k, Mej_solar, Vej0, rho0)


    f = integrate.odeint(dfdt, f0, ate, args)
    #f = rk4(dfdt, f0, ate, args)

    ar = f[:,0]
    au = f[:,1]
    
    
    if not np.isfinite(f).all() or (ar<0).any() or (au<0).any() or (np.diff(ar)/ar[:-1]<-1.0e-10).any() or (np.diff(au)/au[:-1]>1.0e-10).any():
        print("Bad integration")
    """

    ar, au = shock.shockEvolRK4(ate, r0, umax, 
                                Mej_solar*Msun, rho0, Ei, k, umin)


    P = np.zeros(t.shape)

    wopts = None

    for i in range(len(t)):
        amu = c * (ate - t[i]) / ar

        args = (amu, ate, au, ar, nu[i], n0, p, epsE, epsB, ksiN, specType)

        res = integrate.quad(dP, 0.0, np.pi, args, full_output=1, wopts=wopts,
                                epsrel=intrtol)
        P[i] = res[0]
        #print(res)
        #wopts = (res[2]['momcom'], res[2]['chebmo'])

    Fnu = cgs2mJy * math.sqrt(3.0)*ee*ee*ee*(p-1) * P / (
                                                    8*np.pi*dL*dL * me *c*c)

    return Fnu

if __name__ == "__main__":
    
    gmax = 3.5
    umax = np.sqrt(gmax*gmax-1.0)
    umin = 0.1
    q = 5
    Ei = 2.0e51
    Mej = 1.0e-5
    n0 = 8.0e-5
    p = 2.2
    epsE = 1.0e-1
    epsB = 1.0e-2
    ksiN = 1.0
    dL = 1.23e26
    specType = 1

    Y = np.array([umax, umin, Ei, q, Mej, n0, p, epsE, epsB, ksiN, dL])
    t = np.logspace(3, 9, num=100, base=10.0)
    nu = np.empty(t.shape)
    nu[:] = 6.0e9

    Fnu = fluxDensityCocoon(t, nu, 3, specType, *Y)

    import sys
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1,1)
    ax.plot(t, Fnu, 'k-')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t_{obs}$ (s)")
    ax.set_ylabel(r"$F_\nu$ (mJy)")
    
    fig, ax = plt.subplots(1,1)
    ax.plot(t*sec2day, Fnu, 'k-')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1.0, 1.0e3)
    ax.set_xlabel(r"$t_{obs}$ (d)")
    ax.set_ylabel(r"$F_\nu$ (mJy)")

    r0 = 1.0e9
    rho0 = mp * n0
    Vej0 = 4.0/3.0*np.pi*r0*r0*r0

    t0 = r0 / c
    t1 = 1.0e12 * t0
    NT = 15000
    ate = np.logspace(math.log10(t0), math.log10(t1), num=NT, base=10.0)
    f0 = np.array([r0, umax])
    args = (umax, umin, Ei, q, Mej, Vej0, rho0)
    f = integrate.odeint(dfdt, f0, ate, args)
    ar = f[:,0]
    au = f[:,1]

    fig, ax = plt.subplots(2,1)
    ax[0].plot(ate, ar, 'k-')
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel(r"$t_{lab}$ (s)")
    ax[0].set_ylabel(r"$r$ (cm)")
    ax[1].plot(ate, au, 'k-')
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlabel(r"$t_{lab}$ (s)")
    ax[1].set_ylabel(r"$u$")

    intrtol = 1.0e-16
    Fnu0 = fluxDensityCocoon(t, nu, 3, specType, *Y)
    intrtol = 1.0e-14
    Fnu1 = fluxDensityCocoon(t, nu, 3, specType, *Y)
    intrtol = 1.0e-12
    Fnu2 = fluxDensityCocoon(t, nu, 3, specType, *Y)
    intrtol = 1.0e-8
    Fnu3 = fluxDensityCocoon(t, nu, 3, specType, *Y)
    intrtol = 1.0e-6
    Fnu4 = fluxDensityCocoon(t, nu, 3, specType, *Y)
    intrtol = 1.0e-4
    Fnu5 = fluxDensityCocoon(t, nu, 3, specType, *Y)
    intrtol = 1.0e-2
    Fnu6 = fluxDensityCocoon(t, nu, 3, specType, *Y)

    print(np.fabs((Fnu6-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu5-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu4-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu3-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu2-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu1-Fnu0)/Fnu0).max())

    plt.draw()
    plt.show()

    sys.exit()
