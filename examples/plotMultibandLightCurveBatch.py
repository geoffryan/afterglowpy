import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

# Jet Parameters
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

# Time and Frequencies
Nt = 100
Nnu = 3
t = np.empty((Nt, Nnu))
nu = np.empty((Nt, Nnu))

ta = 1.0e-1 * grb.day2sec
tb = 1.0e3 * grb.day2sec
nuR = 6.0e9
nuO = 1.0e14
nuX = 1.0e18

t[:, :] = np.geomspace(ta, tb, num=100)[:, None]
nu[:, 0] = nuR
nu[:, 1] = nuO
nu[:, 2] = nuX

# Calculate!
print("Calculate!")
Fnu = grb.fluxDensity(t, nu, jetType, specType, *Y)

# Plot!
print("Plot")

tday = t * grb.sec2day

fig, ax = plt.subplots(1, 1)
ax.plot(tday[:, 0], Fnu[:, 0], ls='-', label=r'$\nu=6$ GHz')
ax.plot(tday[:, 1], Fnu[:, 1], ls='--', label=r'$\nu=10^{14}$ Hz')
ax.plot(tday[:, 2], Fnu[:, 2], ls='-.', label=r'$\nu=10^{18}$ Hz')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
ax.legend()
fig.tight_layout()

print("Saving lc_multi_batch.png")
fig.savefig("lc_multi_batch.png")
plt.close(fig)
