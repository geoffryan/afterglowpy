import time
import numpy as np
import matplotlib.pyplot as plt
import grbpy as grb
import scipy.integrate as integrate

t = np.logspace(4, 16, 10000)
R0 = t[0] * grb.c
u0 = 10.0

Mej = 1.0e27
rho0 = 1.0e-27
Einj = 1.0e56
k = 5
umin = 1.0

ta = time.clock()
R, u = grb.shockEvolRK4(t, R0, u0, Mej, rho0, Einj, k, umin)
tb = time.clock()
print("CRK4: {0} s".format(tb-ta))

f0 = np.array([R0, u0])
Vej0 = 4./3. * np.pi*R0*R0*R0
args = (u0, umin, Einj, k, Mej/grb.Msun, Vej0, rho0)
ta = time.clock()
foi = integrate.odeint(grb.dfdt, f0, t, args)
tb = time.clock()
print("Odeint: {0} s".format(tb-ta))

Roi = foi[:,0]
uoi = foi[:,1]
ta = time.clock()
frk = grb.rk4(grb.dfdt, f0, t, args)
tb = time.clock()
print("PyRK4: {0} s".format(tb-ta))
Rrk = frk[:,0]
urk = frk[:,1]

fig, ax = plt.subplots(2,2)
ax[0,0].plot(t, R, label="C - RK4")
ax[0,0].plot(t, Roi, label="SciPy odeint")
ax[0,0].plot(t, Rrk, label="Py - RK4")
ax[0,1].plot(t, u)
ax[0,1].plot(t, uoi)
ax[0,1].plot(t, urk)
ax[1,0].plot(t, (Rrk-R)/R)
ax[1,0].plot(t, (Roi-R)/R)
ax[1,1].plot(t, (urk-u)/u)
ax[1,1].plot(t, (uoi-u)/u)

ax[0,0].set_xscale('log')
ax[0,0].set_yscale('log')
ax[0,0].legend()
ax[0,1].set_xscale('log')
ax[0,1].set_yscale('log')
ax[1,0].set_xscale('log')
ax[1,1].set_xscale('log')

print("PyRK4-CRK4   max diff: " + str(np.fabs(((Rrk-R)/R)).max()))
print("odeint-CRK4  max diff: " + str(np.fabs(((Roi-R)/R)).max()))
print("odeint-PyRK4 max diff: " + str(np.fabs(((Roi-Rrk)/Rrk)).max()))
print("PyRK4-CRK4   max diff: " + str(np.fabs(((urk-u)/u)).max()))
print("odeint-CRK4  max diff: " + str(np.fabs(((uoi-u)/u)).max()))
print("odeint-PyRK4 max diff: " + str(np.fabs(((uoi-urk)/urk)).max()))

plt.show()

