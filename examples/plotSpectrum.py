import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

jetType = 4
specType = 0
thV = 0.4
E0 = 1.0e52
thC = 0.1
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

T = 1.0 * grb.day2sec  # spectrum at 1 day

nua = 1.0e6   # Frequencies in Hz
nub = 1.0e20  # Frequencies in Hz

nu = np.geomspace(nua, nub, num=100)
t = np.empty(nu.shape)
t[:] = T

print("Calculating")
Fnu = grb.fluxDensity(t, nu, jetType, specType, *Y)

print("Writing spec.txt")
f = open("spec.txt", 'w')
f.write("# t " + str(t[0]) + ' (s)\n')
f.write("# jetType " + str(jetType) + " specType " + str(specType)+"\n")
f.write("# " + " ".join([str(y) for y in Y]) + "\n")
for i in range(len(t)):
    f.write("{0:.6e} {1:.6e}\n".format(nu[i], Fnu[i]))
f.close()

print("Plotting")
fig, ax = plt.subplots(1, 1)

ax.plot(nu, Fnu)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\nu$ (Hz)')
ax.set_ylabel(r'$F_\nu$[1 day] (mJy)')

fig.tight_layout()
print("Saving figure spec.png")
fig.savefig("spec.png")
plt.close(fig)
