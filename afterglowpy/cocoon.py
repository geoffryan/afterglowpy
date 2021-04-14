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


def dP(costheta, amu, ate, au, ar, nu, n0, p, epsE, epsB, ksiN, specType):

    mu = costheta
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

    em = jet.emissivity(nu, r, mu, te, u, us, n0, p, epsE,
                        epsB, ksiN, specType)

    return 2*np.pi * em


def fluxDensity(t, nu, **kwargs):

    t = np.array(t)

    specType = kwargs['specType']
    uMax = kwargs['uMax']
    uMin = kwargs['uMin']
    Er = kwargs['Er']
    k = kwargs['k']
    MFast_solar = kwargs['MFast_solar']
    n0 = kwargs['n0']
    p = kwargs['p']
    epsilon_e = kwargs['epsilon_e']
    epsilon_B = kwargs['epsilon_B']
    ksiN = kwargs['xi_N']
    dL = kwargs['d_L']

    # Energy injection variables (off by default)
    L0 = kwargs['L0'] if 'L0' in kwargs else 0.0
    q = kwargs['q'] if 'q' in kwargs else 0.0
    ts = kwargs['ts'] if 'ts' in kwargs else 0.0
   
    # Numerical integration variables
    rtol = kwargs['rtol'] if 'rtol' in kwargs else 1.0e-3
    tRes = kwargs['tRes'] if 'tRes' in kwargs else 1000
    latRes = kwargs['latRes'] if 'latRes' in kwargs else 0

    rho0 = mp * n0
    Mej = MFast_solar * Msun
    u0 = uMax
    g0 = math.sqrt(1+u0*u0)
    bes0 = 4*u0*g0 / (4*u0*u0+3)
    Rd = math.pow(9*g0*g0*Mej / (4*np.pi*(g0+1)*(4*u0*u0+3)*rho0), 1./3.)
    td = Rd / (bes0 * c)

    t0 = min(1.0e-2*td, 5.0e-1 * g0*g0*t.min(), 5.0e-1 * t.min()/(1+bes0))
    t1 = 2. * g0*g0*t.max()

    NT = int(tRes * math.log10(t1/t0))
    # print("{0:.3e} {1:.3e} {2:.3e} {3:d}".format(t0, t1, t1/t0, NT))

    r0 = bes0*c*t0

    # Vej0 = 4.0/3.0*np.pi*r0*r0*r0

    ate = np.logspace(math.log10(t0), math.log10(t1), num=NT, base=10.0)

    ar, au = shock.shockEvolRK4(ate, r0, uMax,
                                MFast_solar*Msun, rho0, Er, k, uMin, L0, q, ts)

    P = np.zeros(t.shape)

    wopts = None

    for i in range(len(t)):
        amu = c * (ate - t[i]) / ar

        args = (amu, ate, au, ar, nu[i], n0, p, epsilon_e, epsilon_B, ksiN,
                specType)

        res = integrate.quad(dP, 0.0, 1.0, args, full_output=1, wopts=wopts,
                             epsrel=rtol)
        P[i] = res[0]

    Fnu = cgs2mJy * P / (4*np.pi*dL*dL)

    return Fnu
