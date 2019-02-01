# import sys
# import math
import numpy as np
import matplotlib.pyplot as plt
import grbpy as grb

day = 86400.0

jetType = 4
specType = 0
thV = 0.5
E0 = 1.0e52
thC = 0.01
thW = 0.3
b = 4
L0 = 0.0  # 1.0e47
q = 0.0  # 1.0
ts = 0.0  # 1.0e5
n0 = 1.0e-3
p = 2.15
epse = 1.0e-1
epsB = 1.0e-3
ksiN = 1.0
dL = 1.23e26

Y = np.array([thV, E0, thC, thW, b, L0, q, ts, n0, p, epse, epsB, ksiN, dL])

ta = 1.0e0 * day
tb = 1.0e3 * day
t = np.logspace(np.log10(ta), np.log10(tb), base=10.0, num=100)

nu = np.empty(t.shape)
nu[:] = 1.0e14


print("Calculating")
Fnu = grb.fluxDensity(t, nu, jetType, specType, *Y)

print("Writing lc.txt")
f = open("lc.txt", 'w')
f.write("# nu " + str(nu[0]) + '\n')
f.write("# jetType " + str(jetType) + " specType " + str(specType)+"\n")
f.write("# " + " ".join([str(y) for y in Y]) + "\n")
for i in range(len(t)):
    f.write("{0:.6e} {1:.6e}\n".format(t[i], Fnu[i]))
f.close()


print("Plotting")
title = r"$E_0$ = {0:.01f}x$10^{{52}}$erg".format(E0/1.0e52)
title += r"   $n_0$ = {0:.01f}cm$^{{-3}}$".format(n0)
title += r"   $\nu$ = $10^{{14}}$Hz".format(n0)

fig, ax = plt.subplots(1, 1)
"""
thVs = [0.1, 0.5, 1.0, 1.4]
for thV in thVs:
    Y[0] = thV
    Fnu = grb.fluxDensity(t, nu, jetType, specType, *Y)
    ax.plot(t/day, Fnu, label=r'$\theta_V$ = {0:.02f} rad'.format(thV))
ax.legend()
"""
ax.plot(t/day, Fnu)
Fnu = grb.fluxDensity(t, nu, jetType, specType, *Y, spread=False)
ax.plot(t/day, Fnu)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
ax.set_title(title)
fig.tight_layout()
# ax.set_ylim(1.0e-8, 1.0e-4)
print("Saving figure lc.pdf")
fig.savefig("lc.pdf")
plt.close(fig)
