import time
import numpy as np
import matplotlib.pyplot as plt
import grbpy as grb

NT = 10
NNU = 10
Ntrials = 10
jetTypes = [-1,0,1,2,3,4]

Yj = np.array([0.5, 1.0e52, 0.05, 0.5, 1.0e-3, 2.2, 0.1, 0.01, 1.0, 1.0e26])
Yc = np.array([10.0, 1.0, 1.0e55, 5, 1.0e-7, 1.0e-3, 2.2, 0.1, 0.01, 1.0, 1.0e26])

Ts = np.logspace(0, 3, num=NT, base=10.0) * 86400.0
NUs = np.logspace(3, 20, num=NNU, base=10.0)

t = np.empty(NT * NNU)
nu = np.empty(NT * NNU)

for i in range(NNU):
    t[i*NT:(i+1)*NT] = Ts
    nu[i*NT:(i+1)*NT] = NUs[i]


times = np.empty((len(jetTypes), Ntrials))
Fnu = np.empty((len(jetTypes), NT*NNU))

for i in range(Ntrials):
    print("Trial {0:03d} of {1:03d}".format(i+1, Ntrials))
    for j, jt in enumerate(jetTypes):
        if jt == 3:
            Y = Yc
        else:
            Y = Yj
    
        ta = time.clock()
        Fnu[j] = grb.fluxDensity(t, nu, jt, 0, *Y)
        tb = time.clock()
        times[j,i] = tb-ta

print("Number of evaluation points: {0:d}".format(NT*NNU))
for j, jt in enumerate(jetTypes):
    print("Jet Type {0: d} avg time: {1:.3e} s   ({2:.3e} s per point)".format(jt, times[j].mean(), times[j].mean()/(NT*NNU)))

print("Plotting time_trial_fluxes.png")
fig, ax = plt.subplots(1,1)
for j in range(len(jetTypes)):
    for i in range(NNU):
        ax.plot(Ts, Fnu[j,i*NT:(i+1)*NT])
ax.set_xscale('log')
ax.set_yscale('log')
fig.savefig("time_trial_fluxes.png")

