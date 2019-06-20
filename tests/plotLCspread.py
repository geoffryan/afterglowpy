import sys
# import math
import numpy as np
import matplotlib.pyplot as plt
import grbpy as grb


def f_struct(th, jetType, Y):
    thC = Y[2]
    thW = Y[3]
    b = Y[4]
    if jetType == -2:
        if th > thC and th < thW:
            return 1.0
    elif jetType == -1:
        if th < thC:
            return 1.0
    elif jetType == 0:
        if th < thW:
            return np.exp(-0.5*th*th/(thC*thC))
    else:
        if th < thW:
            return np.power(1.0 + th*th/(thC*thC), -0.5*b)

    return 0.0


day = 86400.0

jetType = 4
specType = 0
thV = 0.5
E0 = 1.0e52
thC = 0.08
thW = 0.4
b = 6
L0 = 0.0  # 1.0e47
q = 0.0  # 1.0
ts = 0.0  # 1.0e5
n0 = 1.0e-3
p = 2.15
epse = 1.0e-1
epsB = 1.0e-3
ksiN = 1.0
dL = 1.23e26
g0 = 0.0

Y = np.array([thV, E0, thC, thW, b, L0, q, ts, n0, p, epse, epsB, ksiN, dL,
              g0])

ta = 1.0e-2 * day
tb = 1.0e4 * day
t = np.logspace(np.log10(ta), np.log10(tb), base=10.0, num=1000)

nu = np.empty(t.shape)
nu[:] = 1.0e14

latRes = 5

Nth = int(latRes * thW / thC)
th = np.linspace(0.0, thW, Nth+1)

Fnu0b = np.empty((Nth, len(t)))
Fnu1b = np.empty((Nth, len(t)))
Fnu2b = np.empty((Nth, len(t)))
Fnu3b = np.empty((Nth, len(t)))
Fnu5b = np.empty((Nth, len(t)))
Fnu7b = np.empty((Nth, len(t)))

print("Calculating")
print("Full spread 0")
sys.stdout.flush()
Fnu0 = grb.fluxDensity(t, nu, jetType, specType, *Y, spread=0, latRes=5)
print("Full spread 1")
sys.stdout.flush()
Fnu1 = grb.fluxDensity(t, nu, jetType, specType, *Y, spread=1, latRes=5)
print("Full spread 2")
sys.stdout.flush()
Fnu2 = grb.fluxDensity(t, nu, jetType, specType, *Y, spread=2, latRes=5)
print("Full spread 4")
sys.stdout.flush()
Fnu3 = grb.fluxDensity(t, nu, jetType, specType, *Y, spread=4, latRes=5)
print("Full spread 5")
sys.stdout.flush()
Fnu5 = grb.fluxDensity(t, nu, jetType, specType, *Y, spread=5, latRes=5)
print("Full spread 7")
sys.stdout.flush()
Fnu7 = grb.fluxDensity(t, nu, jetType, specType, *Y, spread=7, latRes=5)

Yc = Y.copy()
for i in range(Nth):
    tha = th[i]
    thb = th[i+1]
    thth = 0.5*(tha+thb)
    Eth = Y[1] * f_struct(thth, jetType, Y)
    Yc[1] = Eth
    Yc[2] = tha
    Yc[3] = thb
    print("Cone spread 0")
    Fnu0b[i, :] = grb.fluxDensity(t, nu, -2, specType, *Yc,
                                  thetaCoreGlobal=thC, spread=0)
    print("Cone spread 3")
    Fnu1b[i, :] = grb.fluxDensity(t, nu, -2, specType, *Yc,
                                  thetaCoreGlobal=thC, spread=3)
    print("Cone spread 2")
    Fnu2b[i, :] = grb.fluxDensity(t, nu, -2, specType, *Yc,
                                  thetaCoreGlobal=thC, spread=2)
    print("Cone spread 4")
    Fnu3b[i, :] = grb.fluxDensity(t, nu, -2, specType, *Yc,
                                  thetaCoreGlobal=thC, spread=4)
    print("Cone spread 6")
    Fnu5b[i, :] = grb.fluxDensity(t, nu, -2, specType, *Yc,
                                  thetaCoreGlobal=thC, spread=6)
    print("Cone spread 8")
    Fnu7b[i, :] = grb.fluxDensity(t, nu, -2, specType, *Yc,
                                  thetaCoreGlobal=thC, spread=8)

fig, ax = plt.subplots(1, 1, figsize=(12, 9))

for i in range(Nth):
    ax.plot(t/day, Fnu0b[i], color='C0', alpha=0.5)
    ax.plot(t/day, Fnu1b[i], color='C1', alpha=0.5)
    # ax.plot(t/day, Fnu2b[i], color='C2', alpha=0.5)
    # ax.plot(t/day, Fnu3b[i], color='C3', alpha=0.5)
    ax.plot(t/day, Fnu5b[i], color='C4', alpha=0.5)
    ax.plot(t/day, Fnu7b[i], color='C5', alpha=0.5)
ax.plot(t/day, Fnu0b.sum(axis=0), color='C0', ls='--')
ax.plot(t/day, Fnu1b.sum(axis=0), color='C1', ls='--')
ax.plot(t/day, Fnu2b.sum(axis=0), color='C2', ls='--')
ax.plot(t/day, Fnu3b.sum(axis=0), color='C3', ls='--')
ax.plot(t/day, Fnu5b.sum(axis=0), color='C4', ls='--')
ax.plot(t/day, Fnu7b.sum(axis=0), color='C5', ls='--')
ax.plot(t/day, Fnu0, color='C0')
ax.plot(t/day, Fnu1, color='C1')
ax.plot(t/day, Fnu2, color='C2')
ax.plot(t/day, Fnu3, color='C3')
ax.plot(t/day, Fnu5, color='C4')
ax.plot(t/day, Fnu7, color='C5')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
fig.tight_layout()
ax.set_ylim(0.3 * min(Fnu0.min(), Fnu1.min(), Fnu2.min(), Fnu3.min()),
            3.0 * max(Fnu0.max(), Fnu1.max(), Fnu2.max(), Fnu3.max()))
print("Saving figure lc.pdf")
fig.savefig("lc.pdf")
plt.close(fig)

fig, ax = plt.subplots(1, 1)
ax.plot(t/day, Fnu5, color='k')
for i in range(Nth):
    tha = th[i]
    thb = th[i+1]
    thth = 0.5*(tha+thb)
    Eth = Y[1] * f_struct(thth, jetType, Y)
    Yc[1] = Eth
    Yc[2] = tha
    Yc[3] = thb
    Fnu = grb.fluxDensity(t, nu, -2, specType, *Yc,
                          thetaCoreGlobal=thC, spread=6)
    ax.plot(t/day, Fnu, color='C0')
    Yc[2] = thb
    Yc[3] = thb
    Fnu = grb.fluxDensity(t, nu, -1, specType, *Yc,
                          thetaCoreGlobal=thC, spread=6)
    ax.plot(t/day, Fnu, color='C1')
    Fnu = grb.fluxDensity(t, nu, -1, specType, *Yc,
                          thetaCoreGlobal=thC, spread=0)
    ax.plot(t/day, Fnu, color='C2')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
fig.tight_layout()
ax.set_ylim(0.3 * min(Fnu0.min(), Fnu1.min(), Fnu2.min(), Fnu3.min()),
            3.0 * max(Fnu0.max(), Fnu1.max(), Fnu2.max(), Fnu3.max()))
print("Saving figure lcSpreadSection.pdf")
fig.savefig("lcSpreadSection.pdf")
plt.close(fig)
