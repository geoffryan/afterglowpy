import math
import numpy as np
import scipy.integrate as integrate
from . import shock
from . import jet

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
    u = au[ia]*math.pow(te/ate[ia], math.log(au[ib]/au[ia]) 
                                        / math.log(ate[ib]/ate[ia]))
    r = ar[ia]*math.pow(te/ate[ia], math.log(ar[ib]/ar[ia])
                                        / math.log(ate[ib]/ate[ia]))

    g = math.sqrt(u*u+1)

    us = 4*u*g / math.sqrt(8*u*u+9)

    em = jet.emissivity(nu, r, math.sin(theta), mu, te, u, us, n0, p, epsE,
                        epsB, ksiN, specType)

    return 2*np.pi * em


def fluxDensity(t, nu, jetType, specType, umax, umin, Ei, k, Mej_solar, L0, q,
                        n0, p, epsE, epsB, ksiN, dL, tRes=1000, latRes=0,
                        rtol=1.0e-3):

    rho0 = mp * n0
    Mej = Mej_solar * Msun
    u0 = umax
    g0 = math.sqrt(1+u0*u0)
    bes0 = 4*u0*g0 / (4*u0*u0+3)
    Rd = math.pow(9*g0*g0*Mej / (4*np.pi*(g0+1)*(4*u0*u0+3)*rho0), 1./3.)
    td = Rd / (bes0 * c)

    t0 = min(1.0e-2*td, 5.0e-1 * g0*g0*t.min(), 5.0e-1 * t.min()/(1+bes0))
    t1 = 2. * g0*g0*t.max()

    NT = int(tRes * math.log10(t1/t0))
    #print("{0:.3e} {1:.3e} {2:.3e} {3:d}".format(t0, t1, t1/t0, NT))

    r0 = bes0*c*t0

    Vej0 = 4.0/3.0*np.pi*r0*r0*r0

    ate = np.logspace(math.log10(t0), math.log10(t1), num=NT, base=10.0)


    ar, au = shock.shockEvolRK4(ate, r0, umax, 
                                Mej_solar*Msun, rho0, Ei, k, umin, L0, q)


    P = np.zeros(t.shape)

    wopts = None

    for i in range(len(t)):
        amu = c * (ate - t[i]) / ar

        args = (amu, ate, au, ar, nu[i], n0, p, epsE, epsB, ksiN, specType)

        res = integrate.quad(dP, 0.0, np.pi, args, full_output=1, wopts=wopts,
                                epsrel=rtol)
        P[i] = res[0]

    Fnu = cgs2mJy * P / (4*np.pi*dL*dL)

    return Fnu


if __name__ == "__main__":
    
    gmax = 3.5
    umax = np.sqrt(gmax*gmax-1.0)
    umin = 0.1
    k = 5
    Ei = 2.0e51
    Mej = 1.0e-5
    L0 = 0.0
    q = 0.5
    n0 = 8.0e-5
    p = 2.2
    epsE = 1.0e-1
    epsB = 1.0e-2
    ksiN = 1.0
    dL = 1.23e26
    specType = 0

    Y = np.array([umax, umin, Ei, k, Mej, L0, q, n0, p, epsE, epsB, ksiN, dL])
    t = np.logspace(3, 9, num=100, base=10.0)
    nu = np.empty(t.shape)
    nu[:] = 6.0e9

    Fnu = fluxDensity(t, nu, 3, specType, *Y)

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

    t0 = r0 / c
    t1 = 1.0e12 * t0
    NT = 15000
    ate = np.logspace(math.log10(t0), math.log10(t1), num=NT, base=10.0)
    ar, au = shock.shockEvolRK4(ate, r0, umax, 
                                Mej_solar*Msun, rho0, Ei, k, umin, L0, q)

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
    Fnu0 = fluxDensity(t, nu, 3, specType, *Y)
    intrtol = 1.0e-14
    Fnu1 = fluxDensity(t, nu, 3, specType, *Y)
    intrtol = 1.0e-12
    Fnu2 = fluxDensity(t, nu, 3, specType, *Y)
    intrtol = 1.0e-8
    Fnu3 = fluxDensity(t, nu, 3, specType, *Y)
    intrtol = 1.0e-6
    Fnu4 = fluxDensity(t, nu, 3, specType, *Y)
    intrtol = 1.0e-4
    Fnu5 = fluxDensity(t, nu, 3, specType, *Y)
    intrtol = 1.0e-2
    Fnu6 = fluxDensity(t, nu, 3, specType, *Y)

    print(np.fabs((Fnu6-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu5-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu4-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu3-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu2-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu1-Fnu0)/Fnu0).max())

    plt.show()

    sys.exit()
