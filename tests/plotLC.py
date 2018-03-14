import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import grbpy as grb

day = 86400.0

jetType = 0
specType = 0
thV = 0.5
E0 = 1.0e53
thC = 0.08
thW = 0.3
n0 = 1.0e-2
p = 2.15
epse = 1.0e-1
epsB = 1.0e-3
ksiN = 1.0
dL = 1.23e26

Y = np.array([thV, E0, thC, thW, n0, p, epse, epsB, ksiN, dL])

ta = 1.0e0 * day
tb = 1.0e3 * day
t = np.logspace(np.log10(ta), np.log10(tb), base=10.0, num=100)

nu  = np.empty(t.shape)
nu[:] = 1.0e18


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
fig, ax = plt.subplots(1,1)
ax.plot(t/day, Fnu)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
print("Saving figure lc.png")
fig.savefig("lc.png")
plt.close(fig)
