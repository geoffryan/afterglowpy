import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import grbpy as grb
import grbpy.shock as shock


if __name__ == "__main__":
    
    gmax = 3.5
    umax = np.sqrt(gmax*gmax-1.0)
    umin = 0.1
    k = 5
    Ei = 2.0e51
    Mej_solar = 1.0e-5
    L0 = 1.0e45
    q = 0.5
    n0 = 8.0e-5
    p = 2.2
    epsE = 1.0e-1
    epsB = 1.0e-2
    ksiN = 1.0
    dL = 1.23e26
    specType = 0

    Y = np.array([umax, umin, Ei, k, Mej_solar, L0, q, n0, p, epsE, epsB, 
                    ksiN, dL])
    t = np.logspace(3, 9, num=100, base=10.0)
    nu = np.empty(t.shape)
    nu[:] = 6.0e9

    Fnu = grb.fluxDensity(t, nu, 3, specType, *Y)

    fig, ax = plt.subplots(1,1)
    ax.plot(t, Fnu, 'k-')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t_{obs}$ (s)")
    ax.set_ylabel(r"$F_\nu$ (mJy)")
    
    fig, ax = plt.subplots(1,1)
    ax.plot(t*grb.sec2day, Fnu, 'k-')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1.0, 1.0e3)
    ax.set_xlabel(r"$t_{obs}$ (d)")
    ax.set_ylabel(r"$F_\nu$ (mJy)")

    r0 = 1.0e9
    rho0 = grb.mp * n0

    t0 = r0 / grb.c
    t1 = 1.0e12 * t0
    NT = 15000
    ate = np.logspace(math.log10(t0), math.log10(t1), num=NT, base=10.0)
    ar, au = shock.shockEvolRK4(ate, r0, umax, 
                                Mej_solar*grb.Msun, rho0, Ei, k, umin, L0, q)

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
    Fnu0 = grb.fluxDensity(t, nu, 3, specType, *Y)
    intrtol = 1.0e-14
    Fnu1 = grb.fluxDensity(t, nu, 3, specType, *Y)
    intrtol = 1.0e-12
    Fnu2 = grb.fluxDensity(t, nu, 3, specType, *Y)
    intrtol = 1.0e-8
    Fnu3 = grb.fluxDensity(t, nu, 3, specType, *Y)
    intrtol = 1.0e-6
    Fnu4 = grb.fluxDensity(t, nu, 3, specType, *Y)
    intrtol = 1.0e-4
    Fnu5 = grb.fluxDensity(t, nu, 3, specType, *Y)
    intrtol = 1.0e-2
    Fnu6 = grb.fluxDensity(t, nu, 3, specType, *Y)

    print(np.fabs((Fnu6-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu5-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu4-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu3-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu2-Fnu0)/Fnu0).max())
    print(np.fabs((Fnu1-Fnu0)/Fnu0).max())

    plt.show()

    sys.exit()
