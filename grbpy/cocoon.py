import math
import numpy as np
import scipy.integrate as integrate

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


def dfdt(f, t, Ei, ts, q, Mej, Vej0, rho0):

    R = f[0]
    u = f[1]

    g = math.sqrt(u*u+1)

    dRdt = 4*c* u * g / (4*u*u + 3)

    dEdt = 0.0
    if ts > 0.0 and t < ts:
        dEdt = Ei/ts * math.pow(t/ts, q)

    A_g1 = 4.0/3.0 * np.pi * rho0 * R*R/(g*g) * (4*g*g*g*g - 5*g*g + 1) * dRdt
    A_g2 = (Mej - 4.0/3.0 * rho0 * Vej0 * (19*g*g + 4*g)
            + 8.0/9.0*np.pi*rho0 * R*R*R * (4*g*g*g*g - 5*g*g*g + 5*g*g - 1.0)
                / (g*g*g))
    dgdt = (dEdt/(c*c) - A_g1) / A_g2

    dudt = g*dgdt/u

    return np.array([dRdt,dudt])

def dP(theta, amu, ate, aus, ar, nu, n0, p, epsE, epsB, ksiN):

    mu = math.cos(theta)
    ib = np.searchsorted(amu, mu)
    ia = ib-1

    N = amu.shape[0]
   
    te = ((mu-amu[ia])*ate[ib] + (amu[ib]-mu)*ate[ia]) / (amu[ib]-amu[ia])
    u = ((mu-amu[ia])*aus[ib] + (amu[ib]-mu)*aus[ia]) / (amu[ib]-amu[ia])
    r = ((mu-amu[ia])*ar[ib] + (amu[ib]-mu)*ar[ia]) / (amu[ib]-amu[ia])

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
    gc = 6*np.pi * g * me*c / (sigmaT * B*B * te)
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

def fluxDensityCocoon(t, nu, umax, Ei, ts, q, Mej, n0, p, epsE, epsB, 
                        ksiN, dL):

    r0 = 1.0e9
    rho0 = mp * n0

    Vej0 = 4.0/3.0*np.pi*r0*r0*r0

    #t0 = 1.0e-2 * day2sec
    t0 = r0 / c
    #t1 = 1.0e7 * day2sec
    t1 = 1.0e10 * t0
    NT = 12000
    ate = np.logspace(math.log10(t0), math.log10(t1), num=NT, base=10.0)

    f0 = np.array([r0, umax])
    args = (Ei, ts, q, Mej, Vej0, rho0)

    f = integrate.odeint(dfdt, f0, ate, args)

    ar = f[:,0]
    au = f[:,1]

    P = np.zeros(t.shape)

    wopts = None

    for i in range(len(t)):
        amu = c * (ate - t[i]) / ar

        args = (amu, ate, au, ar, nu[i], n0, p, epsE, epsB, ksiN)

        res = integrate.quad(dP, 0.0, np.pi, args, full_output=1, wopts=wopts)
        P[i] = res[0]
        #print(res)
        #wopts = (res[2]['momcom'], res[2]['chebmo'])

    Fnu = cgs2mJy * math.sqrt(3.0)*ee*ee*ee*(p-1) * P / (
                                                    8*np.pi*dL*dL * me *c*c)

    return Fnu

if __name__ == "__main__":
    
    gmax = 3.5
    umax = np.sqrt(gmax*gmax-1.0)
    Ei = 1.0e49
    ts = 1.0e7
    q = 0.0
    Mej = 1.0e-8 * Msun
    n0 = 1.0e-4
    p = 2.2
    epsE = 1.0e-1
    epsB = 1.0e-2
    ksiN = 1.0
    dL = 1.23e26

    Y = np.array([umax, Ei, ts, q, Mej, n0, p, epsE, epsB, ksiN, dL])
    t = np.logspace(3, 7, num=100, base=10.0)
    nu = np.empty(t.shape)
    nu[:] = 6.0e9

    Fnu = fluxDensityCocoon(t, nu, *Y)

    import sys
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1,1)
    ax.plot(t, Fnu, 'k-')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$F_\nu$ (mJy)")

    plt.draw()
    plt.show()

    sys.exit()
