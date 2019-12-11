import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

jetType = -1
specType = 0
thV = 0.0
E0 = 1.0e52
thC = 0.2
thW = 0.6
b = 6
L0 = 0.0
q = 0.0
ts = 0.0
n0 = 1.0e-3
p = 2.15
epse = 1.0e-1
epsB = 1.0e-2
ksiN = 1.0
dL = 1.23e26

Y = np.array([thV, E0, thC, thW, b, L0, q, ts, n0, p, epse, epsB, ksiN, dL])

ta = 1.0e-3 * grb.day2sec
tb = 1.0e2 * grb.day2sec

t = np.geomspace(ta, tb, num=200)
nu = 2.0e14

print("Calculating")
Fnu = grb.fluxDensity(t, nu, jetType, specType, *Y, spread=False, latRes=5,
                      g0=100)

print("Writing lc.txt")
f = open("lc.txt", 'w')
f.write("# nu " + str(nu) + '\n')
f.write("# jetType " + str(jetType) + " specType " + str(specType)+"\n")
f.write("# " + " ".join([str(y) for y in Y]) + "\n")
for i in range(len(t)):
    f.write("{0:.6e} {1:.6e}\n".format(t[i], Fnu[i]))
f.close()

print("Plotting")
fig, ax = plt.subplots(1, 1)

ax.plot(t/grb.day2sec, Fnu)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$[1keV] (mJy)')

fig.tight_layout()
print("Saving figure lc.png")
fig.savefig("lc.png")
plt.close(fig)
