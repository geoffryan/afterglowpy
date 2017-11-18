import sys
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import scipy.optimize as opt
import corner
import grbpy as grb
import fit
import getOKCData as data

bounds = fit.bounds

file = sys.argv[1]

f = h5.File(file, "r")
chain = f['chain'][...]
lnprobability = f['lnprobability'][...]
X0 = f['X0'][...]
jetType = f['jetType'][0]
fitVars = f['fitVars'][...]
labels = f['labels'][...]
T = f['t'][...]
NU = f['nu'][...]
FNU = f['Fnu'][...]
FERR = f['eFnu'][...]
UL = f['ul'][...]
INST = f['inst'][...]
f.close()

nwalkers = chain.shape[0]
nsteps = chain.shape[1]
ndim = chain.shape[2]

flatlnprobability = lnprobability.reshape((-1,))
flatchain = chain.reshape((-1,ndim))

imax = np.argmax(flatlnprobability)
x1 = flatchain[imax]

cutChain = chain[:,int(0.6*nsteps):,:]
cutSamples = cutChain[:,:,:].reshape((-1,ndim))
fig = corner.corner(cutSamples, labels=labels[fitVars])
fig.savefig("corner_noburn.png")
plt.close(fig)

N = len(T)

X1 = X0.copy()
X1[fitVars] = x1[:]

lpargs=(fit.logPriorFlat, fit.logLikeChi2, jetType, 
            X1, fitVars, None, T, NU, FNU, FERR, 
            False)
lpargsOpt=(fit.logPriorFlat, fit.logLikeChi2, jetType, 
            X1, fitVars, None, T, NU, FNU, FERR, 
            True)

f = h5.File(file, "r")
if "minChi2_val" in f.keys() and "minChi2_x" in f.keys():
    x2 = f['minChi2_x'][...]
    f.close()
else:
    f.close()
    fit.printLP = True
    res = opt.minimize(fit.logpost, x1, args=lpargsOpt, bounds=bounds[fitVars])
    fit.printLP = False
    x2 = res.x[:]

    f = h5.File(file, "a")
    f.create_dataset("minChi2_val", data=res.fun)
    f.create_dataset("minChi2_x", data=res.x)
    f.close()

X2 = X0.copy()
X2[fitVars] = x2

lp1 = flatlnprobability[imax]
redChi21 = -2*lp1 / (ndim-1.0)
print("best MCMC chi2/dof = " + str(redChi21))

lp2 = fit.logpost(x2, *lpargs)
redChi22 = -2*lp2 / (ndim-1.0)
print("best OPT chi2/dof = " + str(redChi22))

fig, ax = plt.subplots(1,1)

t0 = 1*fit.day
t1 = 300*fit.day

t = np.logspace(np.log10(t0), np.log10(t1), num=200, base=10.0)
nu = np.empty(t.shape)

nuX = NU[NU>1.0e16].mean()

X = X0.copy()

for i in range(nwalkers):
    X[fitVars] = chain[i,-1,:]
    Y = fit.getEvalForm(X)
    nu[:] = 6.0e9
    Fnu = grb.fluxDensity(t, nu, jetType, *Y)
    ax.plot(t/fit.day, Fnu, color='k', ls='-', alpha=0.2)
    nu[:] = nuX
    Fnu = grb.fluxDensity(t, nu, jetType, *Y)
    ax.plot(t/fit.day, Fnu, color='k', ls='-', alpha=0.2)

Y1 = fit.getEvalForm(X1)
nu[:] = 6.0e9
Fnu = grb.fluxDensity(t, nu, jetType, *Y1)
ax.plot(t/fit.day, Fnu, color='k', ls='-')
nu[:] = nuX.mean()
Fnu = grb.fluxDensity(t, nu, jetType, *Y1)
ax.plot(t/fit.day, Fnu, color='k', ls='-')

Y2 = fit.getEvalForm(X2)
nu[:] = 6.0e9
Fnu = grb.fluxDensity(t, nu, jetType, *Y2)
ax.plot(t/fit.day, Fnu, color='r', ls='-')
nu[:] = nuX.mean()
Fnu = grb.fluxDensity(t, nu, jetType, *Y2)
ax.plot(t/fit.day, Fnu, color='r', ls='-')

ax.set_xlim(t0/fit.day, t1/fit.day)
ax.set_ylim(1.0e-9, 1.0e-1)

fit.plot_data(ax, T, NU, FNU, FERR, UL, INST)
fig.savefig("lc_all.png")

