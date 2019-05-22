# import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import grbpy as grb


def mu_th(th, thobs, phi):
    return np.cos(th)*np.cos(thobs) + np.sin(th)*np.sin(thobs)*np.cos(phi)


t0 = 1.0e3
t1 = 1.0e12
E0 = 1.0e53
n0 = 1.0
rho0 = n0*grb.mp
tN = np.power(9*E0/(16*np.pi*rho0*grb.c**5), 1.0/3.0)
u0 = np.power(t0/tN, -1.5)
th0 = 0.01
R0 = t0*grb.c * (1-1.0 / (4*u0)) * (1 + 1.0/(4*u0))

# Approx jet break time
tb = 0.4*tN * np.power(th0, 8.0/3.0)

t = np.geomspace(t0, t1, num=10000)

R, u, thj = grb.shock.shockEvolSpreadRK4(t, R0, u0, th0, 0.0, rho0, 0.0, 0.0,
                                         0.0, 0.0, 0.0, 0.0, True)
RN, uN, thjN = grb.shock.shockEvolSpreadRK4(t, R0, u0, th0, 0.0, rho0, 0.0,
                                            0.0, 0.0, 0.0, 0.0, 0.0, False)
tobs_OA = t - R/grb.c
tobsN_OA = t - RN/grb.c

NU = 1.0e+30
p = 2.2
epse = 0.1
epsB = 0.01
xiN = 1.0
dL = 1.0e26
Y = np.array([0.0, E0, th0, th0, 0.0, 0.0, 0.0, 0.0, n0, p, epse, epsB, xiN,
              dL])
nufac = 1.1
NUa = NU/nufac
NUb = NU*nufac

us = 4*u*np.sqrt((u*u+1)/(8*u*u+9))
usN = 4*uN*np.sqrt((uN*uN+1)/(8*uN*uN+9))

enu = np.array([grb.jet.emissivity(NU, R[i], 1.0, 1.0, t[i], u[i], us[i],
                n0, p, epse, epsB, xiN, 0)
                for i in range(len(t))])
enuN = np.array([grb.jet.emissivity(NU, RN[i], 1.0, 1.0, t[i], uN[i], usN[i],
                 n0, p, epse, epsB, xiN, 0)
                 for i in range(len(t))])
enuna = np.array([grb.jet.emissivity(NUa, R[i], 1.0, 1.0, t[i], u[i],
                  us[i], n0, p, epse, epsB, xiN, 0)
                  for i in range(len(t))])
enunb = np.array([grb.jet.emissivity(NUb, R[i], 1.0, 1.0, t[i], u[i],
                  us[i], n0, p, epse, epsB, xiN, 0)
                  for i in range(len(t))])
enunaN = np.array([grb.jet.emissivity(NUa, RN[i], 1.0, 1.0, t[i], uN[i],
                   usN[i], n0, p, epse, epsB, xiN, 0)
                   for i in range(len(t))])
enunbN = np.array([grb.jet.emissivity(NUb, RN[i], 1.0, 1.0, t[i], uN[i],
                   usN[i], n0, p, epse, epsB, xiN, 0)
                   for i in range(len(t))])
beta = np.log(enunb/enuna) / np.log(NUb/NUa)
betaN = np.log(enunbN/enunaN) / np.log(NUb/NUa)

enu *= grb.cgs2mJy / (4*np.pi*dL*dL)
enuN *= grb.cgs2mJy / (4*np.pi*dL*dL)

fig, ax = plt.subplots(3, 2, figsize=(6, 8))
ax[0, 0].plot(t, R)
ax[0, 1].plot(t, u)
ax[1, 0].plot(t, thj)
ax[1, 1].plot(t, tobs_OA)
ax[0, 0].plot(t, RN)
ax[0, 1].plot(t, uN)
ax[1, 0].plot(t, thjN)
ax[1, 1].plot(t, tobsN_OA)
ax[2, 0].plot(t, enu)
ax[2, 0].plot(t, enuN)
ax[2, 1].plot(t, beta)
ax[2, 1].plot(t, betaN)

ax[0, 0].set_xscale('log')
ax[0, 1].set_xscale('log')
ax[1, 0].set_xscale('log')
ax[1, 1].set_xscale('log')
ax[2, 0].set_xscale('log')
ax[2, 1].set_xscale('log')

ax[0, 0].set_yscale('log')
ax[0, 1].set_yscale('log')
# ax[1, 0].set_yscale('log')
ax[1, 1].set_yscale('log')
# ax[1, 1].set_ylim(t0, t1)
ax[2, 0].set_yscale('log')

tobs = np.geomspace(0.01 * tb, 1.0e7, 1000)
nu = np.empty(tobs.shape)
nu[:] = NU
Fnu = grb.fluxDensity(tobs, nu, -1, 0, *Y, spread=True)
FnuN = grb.fluxDensity(tobs, nu, -1, 0, *Y, spread=False)
TH = np.linspace(0.0, 3.0*th0, 40)

enuth = np.empty((len(TH), len(t)))
enuthN = np.empty((len(TH), len(t)))
for i in range(len(TH)):
    for j in range(len(t)):
        enuth[i, j] = grb.jet.emissivity(NU, R[j], 1.0, math.cos(TH[i]),
                                         t[j], u[j], us[j], n0, p, epse, epsB,
                                         xiN, 0)
        if TH[i] > thj[j]:
            enuth[i, j] = 0.0
        enuthN[i, j] = grb.jet.emissivity(NU, RN[j], 1.0, math.cos(TH[i]),
                                          t[j], uN[j], usN[j], n0, p, epse,
                                          epsB, xiN, 0)
        if TH[i] > thjN[j]:
            enuthN[i, j] = 0.0

enuth *= grb.cgs2mJy / (4*np.pi*dL*dL)
enuthN *= grb.cgs2mJy / (4*np.pi*dL*dL)

tobsth = t[None, :] - np.cos(TH)[:, None]*R[None, :] / grb.c
tobsthN = t[None, :] - np.cos(TH)[:, None]*RN[None, :] / grb.c

enuth_i = np.empty(enuth.shape)
enuthN_i = np.empty(enuthN.shape)
for i in range(len(TH)):
    enuth_i[i, :] = np.interp(tobs_OA, tobsth[i], enuth[i])
    enuthN_i[i, :] = np.interp(tobsN_OA, tobsthN[i], enuthN[i])

f = np.sin(TH)[:, None] * enuth_i[:, :]
Fnuth = 2*np.pi * (0.5*(f[:-1, :]+f[1:, :]) * (TH[1:]-TH[:-1])[:, None]
                   ).sum(axis=0)
f = np.sin(TH)[:, None] * enuthN_i[:, :]
FnuthN = 2*np.pi * (0.5*(f[:-1, :]+f[1:, :]) * (TH[1:]-TH[:-1])[:, None]
                    ).sum(axis=0)

# Flux interpolated to t-grid
Fnu_i = np.interp(tobs_OA, tobs, Fnu)
FnuN_i = np.interp(tobsN_OA, tobs, FnuN)
# Effective angular size of patch, w.r.t. emissivity of core
dOm_i = Fnu_i / enu
dOmN_i = FnuN_i / enuN

# One-sided angular size of jet
dOmj = 2*np.pi*(1-np.cos(thj))
dOmjN = 2*np.pi*(1-np.cos(thjN))

# Effective angular size of patch from emissivity integral.
dOmth = Fnuth / enu
dOmthN = FnuthN / enuN

# Approx angular size of beamed patch
th_u = 0.75/np.sqrt(u*u+1)
th_uN = 0.75/np.sqrt(uN*uN+1)
dOm_u = 2*np.pi*(1-np.cos(th_u))
dOm_uN = 2*np.pi*(1-np.cos(th_uN))

iN = len(t) - np.searchsorted(u[::-1], 1.0)
iNN = len(t) - np.searchsorted(uN[::-1], 1.0)
tN_e = tobs_OA[iN]
tNN_e = tobsN_OA[iNN]

fig3, ax3 = plt.subplots(2, 1, figsize=(6, 8))
ax3[0].axvline(tb, ls='-', color='grey')
ax3[0].axvline(tN_e, ls='--', color='grey')
ax3[0].plot(tobs, Fnu, color='C0')
ax3[0].plot(tobs_OA, enu*dOmj, color='C1')
ax3[0].plot(tobs_OA, Fnuth, color='C2')
ax3[0].axvline(tNN_e, ls='--', color='grey', alpha=0.5)
ax3[0].plot(tobs, FnuN, color='C0', alpha=0.5, lw=2)
ax3[0].plot(tobsN_OA, enuN*dOmjN, color='C1', alpha=0.5, lw=2)
ax3[0].plot(tobsN_OA, FnuthN, color='C2', alpha=0.5, lw=2)

ax3[1].axvline(tb, ls='-', color='grey')
ax3[1].axvline(tN_e, ls='--', color='grey')
ax3[1].plot(tobs_OA, dOm_i, color='C0')
ax3[1].plot(tobs_OA, dOmj, color='C1')
ax3[1].plot(tobs_OA, dOmth, color='C2')
ax3[1].axvline(tNN_e, ls='--', color='grey', alpha=0.5)
ax3[1].plot(tobsN_OA, dOmN_i, color='C0', alpha=0.5, lw=2)
ax3[1].plot(tobsN_OA, dOmjN, color='C1', alpha=0.5, lw=2)
ax3[1].plot(tobsN_OA, dOmthN, color='C2', alpha=0.5, lw=2)


ax3[0].set_xlim(tobs[0], tobs[-1])
ax3[0].set_ylim(0.1*Fnu.min(), 10*Fnu.max())
ax3[1].set_xlim(tobs[0], tobs[-1])
ax3[1].set_ylim(th0*th0, 4*np.pi)

ax3[0].set_xscale('log')
ax3[0].set_yscale('log')
ax3[1].set_xscale('log')
ax3[1].set_yscale('log')

ax3[0].set_ylabel(r"$F_{\nu}$ (mJy)")
ax3[1].set_ylabel(r"$\Delta \Omega$ (sr)")
ax3[1].set_xlabel(r"$t_{\mathrm{obs}}$ (s)")

fig4, ax4 = plt.subplots(2, 1, figsize=(6, 8))
TOBS = np.geomspace(0.3*tb, 30*tb, 10)
for i in range(len(TOBS)):
    it = np.searchsorted(tobs_OA, TOBS[i])
    col = 'C{0:d}'.format(i)
    ax4[0].plot(TH, enuth_i[:, it], color=col)
    ax4[0].plot(TH, enuthN_i[:, it], color=col, alpha=0.5)

thje = np.empty(len(tobs))
thjeN = np.empty(len(tobs))
for i in range(len(tobs)):
    thje[i] = grb.jet.find_jet_edge(t, R, thj, tobs[i], 0.0, 0.0, th0, 1.0)
    thjeN[i] = grb.jet.find_jet_edge(t, RN, thjN, tobs[i], 0.0, 0.0, th0, 1.0)

thjenu = np.empty(t.shape)
thjenuN = np.empty(t.shape)
for i in range(len(t)):
    thjenu[i] = TH[enuth_i[:, i] > 0.0][-1]
    thjenuN[i] = TH[enuthN_i[:, i] > 0.0][-1]

ax4[1].plot(tobs, thje, color='C0')
ax4[1].plot(tobs, thjeN, color='C0', alpha=0.5)
ax4[1].plot(tobs_OA, thj, color='C1')
ax4[1].plot(tobsN_OA, thjN, color='C1', alpha=0.5)
ax4[1].plot(tobs_OA, thjenu, color='C2')
ax4[1].plot(tobsN_OA, thjenuN, color='C2', alpha=0.5)

ax4[0].set_yscale('log')
ax4[1].set_xscale('log')
ax4[1].set_yscale('log')
ax4[1].set_xlim(tobs[0], tobs[-1])
ax4[1].set_ylim(None, 0.5*np.pi)

# plt.show()
# sys.exit()

# th = np.linspace(0.0, 0.5*np.pi, 300)
th = np.linspace(0.0, 0.1, 300)

# thobss = [0.0, 0.3, 0.6, 1.0, 0.5*np.pi]
thobss = [0.0]
phi = [0.0, np.pi/3.0, 2*np.pi/3.0, np.pi]
# tobss = [1.0e4, 3.0e5, 1.0e7, 3.0e8, 1.0e10]
tobss = [3.0e1, 1.0e2, 3.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4]

colors = ['C{0:d}'.format(i) for i in range(10)]
lss = ['-', '--', '-.', ':']

for i in range(len(thobss)):
    fig2, ax2 = plt.subplots(1, 3, figsize=(12, 4))
    for j in range(len(phi)):
        ax2[0].plot(th, mu_th(th, thobss[i], phi[j]), color='k', ls=lss[j])
    for j in range(len(tobss))[::-1]:
        muE = grb.c*(t-tobss[j]) / R
        ax2[0].plot(thj, muE, color=colors[j])
        ax2[1].plot(t, muE, color=colors[j])
    ax2[2].plot(th, th, color='k')
    for j in range(len(phi)):
        muth = mu_th(th, thobss[i], phi[j])
        for k in range(len(tobss)):
            muE = grb.c*(t-tobss[k]) / R
            thji = np.interp(muth, muE, thj)
            ax2[2].plot(th, thji, color=colors[k], ls=lss[j])

    ax2[0].set_xlim(th.min(), th.max())
    ax2[0].set_ylim(-1, 1)
    ax2[1].set_ylim(-1, 1)
    ax2[1].set_xscale('log')
    ax2[2].set_xlim(th.min(), th.max())
    ax2[2].set_ylim(th.min(), th.max())

plt.show()
