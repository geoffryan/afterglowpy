import numpy as np
import matplotlib.pyplot as plt
import grbpy as grb

gmax = 10.0
umax = np.sqrt(gmax*gmax-1.0)
umin = 1.0
Ei = 1.0e49
k = 0.0
Mej = 1.0e-8 * grb.Msun
n0 = 1.0e-4
p = 2.2
epsE = 1.0e-1
epsB = 1.0e-2
ksiN = 1.0
dL = 1.23e26

jetType = 3
specType = 0

Y = np.array([umax, umin, Ei, k, Mej, 0.0, 0.0, 0.0, n0, p, epsE, epsB, ksiN,
              dL])
t = np.logspace(3, 7, num=100, base=10.0)
nu = np.empty(t.shape)
nu[:] = 6.0e9

Fnu1 = grb.fluxDensity(t, nu, jetType, specType, *Y)
Fnu2 = grb.fluxDensity(t, nu, jetType, specType, *Y)
Fnu3 = grb.fluxDensity(t, nu, jetType, specType, *Y)
Fnu4 = grb.fluxDensity(t, nu, jetType, specType, *Y)
Fnu5 = grb.fluxDensity(t, nu, jetType, specType, *Y)
Fnu6 = grb.fluxDensity(t, nu, jetType, specType, *Y)
Fnu7 = grb.fluxDensity(t, nu, jetType, specType, *Y)
Fnu8 = grb.fluxDensity(t, nu, jetType, specType, *Y)
Fnu9 = grb.fluxDensity(t, nu, jetType, specType, *Y)
Fnu10 = grb.fluxDensity(t, nu, jetType, specType, *Y)
Fnu11 = grb.fluxDensity(t, nu, jetType, specType, *Y)
Fnu12 = grb.fluxDensity(t, nu, jetType, specType, *Y)

fig, ax = plt.subplots(1, 1)
ax.plot(t, Fnu1, 'b-')
ax.plot(t, Fnu2, 'g-')
ax.plot(t, Fnu3, 'r-')
ax.plot(t, Fnu4, 'c-')
ax.plot(t, Fnu5, 'b--')
ax.plot(t, Fnu6, 'g--')
ax.plot(t, Fnu7, 'r--')
ax.plot(t, Fnu8, 'c--')
ax.plot(t, Fnu9, 'b:')
ax.plot(t, Fnu10, 'g:')
ax.plot(t, Fnu11, 'r:')
ax.plot(t, Fnu12, 'c:')
ax.set_xscale('log')
ax.set_yscale('log')

plt.show()
