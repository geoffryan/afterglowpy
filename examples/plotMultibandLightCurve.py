import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

jetType = 0
specType = 0
thV = 0.5
E0 = 1.0e53
thC = 0.08
thW = 0.3
b = 0
n0 = 1.0e-2
p = 2.15
epse = 1.0e-1
epsB = 1.0e-2
ksiN = 1.0
dL = 1.23e26

Y = np.array([thV, E0, thC, thW, b, 0, 0, 0, n0, p, epse, epsB, ksiN, dL])

ta = 1.0e-1 * grb.day2sec
tb = 1.0e3 * grb.day2sec
t = np.geomspace(ta, tb, num=100)

nuR = 6.0e9
nuO = 1.0e14
nuX = 1.0e18

print("Calc Radio")
FnuR = grb.fluxDensity(t, nuR, jetType, specType, *Y)
print("Calc Optical")
FnuO = grb.fluxDensity(t, nuO, jetType, specType, *Y)
print("Calc X-ray")
FnuX = grb.fluxDensity(t, nuX, jetType, specType, *Y)

print("Plot")

tday = t * grb.sec2day

fig, ax = plt.subplots(1, 1)
ax.plot(tday, FnuR, ls='-', label=r'$\nu=6$ GHz')
ax.plot(tday, FnuO, ls='--', label=r'$\nu=10^{14}$ Hz')
ax.plot(tday, FnuX, ls='-.', label=r'$\nu=10^{18}$ Hz')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
ax.legend()
fig.tight_layout()

print("Saving lc_multi.png")
fig.savefig("lc_multi.png")
plt.close(fig)
