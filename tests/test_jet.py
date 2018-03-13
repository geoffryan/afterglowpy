import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import grbpy as grb


thV = 0.5
E0 = 1.0e52
thC = 0.05
thW = 0.3
n0 = 1.0e-3
p = 2.15
epse = 1.0e-1
epsB = 1.0e-3
ksiN = 1.0
dL = 1.23e26

specType = 0

day = 86400

Y = np.array([thV, E0, thC, thW, n0, p, epse, epsB, ksiN, dL])

ta = 1.0 * day
tb = 1.0e3 * day
t = np.logspace(np.log10(ta), np.log10(tb), base=10.0, num=300)

nuR  = np.empty(t.shape)
nuO  = np.empty(t.shape)
nuX  = np.empty(t.shape)

nuR[:] = 6.0e9
nuO[:] = 1.0e14
nuX[:] = 1.0e18

print("Calc Radio")
FnuRTH = grb.fluxDensity(t, nuR, -1, specType, *Y)
FnuRGA = grb.fluxDensity(t, nuR, 0, specType, *Y)
FnuRPL = grb.fluxDensity(t, nuR, 1, specType, *Y)
FnuRGC = grb.fluxDensity(t, nuR, 2, specType, *Y)
FnuRPS = grb.fluxDensity(t, nuR, 4, specType, *Y)
print("Calc Optical")
FnuOTH = grb.fluxDensity(t, nuO, -1, specType, *Y)
FnuOGA = grb.fluxDensity(t, nuO, 0, specType, *Y)
FnuOPL = grb.fluxDensity(t, nuO, 1, specType, *Y)
FnuOGC = grb.fluxDensity(t, nuO, 2, specType, *Y)
FnuOPS = grb.fluxDensity(t, nuO, 4, specType, *Y)
print("Calc X-ray")
FnuXTH = grb.fluxDensity(t, nuX, -1, specType, *Y)
FnuXGA = grb.fluxDensity(t, nuX, 0, specType, *Y)
FnuXPL = grb.fluxDensity(t, nuX, 1, specType, *Y)
FnuXGC = grb.fluxDensity(t, nuX, 2, specType, *Y)
FnuXPS = grb.fluxDensity(t, nuX, 4, specType, *Y)

c = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
ls = ["-", "--", ":"]

td = t/day

print("Plot Radio")
fig, ax = plt.subplots(1,1)
ax.plot(td, FnuRTH, color=c[0], label="Top Hat")
ax.plot(td, FnuRGA, color=c[1], label="Gaussian")
ax.plot(td, FnuRGC, color=c[2], label="GaussianCore")
ax.plot(td, FnuRPL, color=c[3], label="PowerlawCore")
ax.plot(td, FnuRPS, color=c[4], label="PowerlawSmooth")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
ax.legend()
fig.savefig("jets_radio.png")
plt.close(fig)

print("Plot Optical")
fig, ax = plt.subplots(1,1)
ax.plot(td, FnuOTH, color=c[0], label="Top Hat")
ax.plot(td, FnuOGA, color=c[1], label="Gaussian")
ax.plot(td, FnuOGC, color=c[2], label="GaussianCore")
ax.plot(td, FnuOPL, color=c[3], label="PowerlawCore")
ax.plot(td, FnuOPS, color=c[4], label="PowerlawSmooth")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
ax.legend()
fig.savefig("jets_optical.png")
plt.close(fig)

print("Plot X-ray")
fig, ax = plt.subplots(1,1)
ax.plot(td, FnuXTH, color=c[0], label="Top Hat")
ax.plot(td, FnuXGA, color=c[1], label="Gaussian")
ax.plot(td, FnuXGC, color=c[2], label="GaussianCore")
ax.plot(td, FnuXPL, color=c[3], label="PowerlawCore")
ax.plot(td, FnuXPS, color=c[4], label="PowerlawSmooth")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
ax.legend()
fig.savefig("jets_xray.png")
plt.close(fig)

print("Plot tophat")
fig, ax = plt.subplots(1,1)
ax.plot(td, FnuRTH, color=c[0], ls=ls[0], label="Radio")
ax.plot(td, FnuOTH, color=c[0], ls=ls[1], label="Optical")
ax.plot(td, FnuXTH, color=c[0], ls=ls[2], label="X-Ray")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
ax.legend()
fig.savefig("jets_tophat.png")
plt.close(fig)

print("Plot Gaussian")
fig, ax = plt.subplots(1,1)
ax.plot(td, FnuRGA, color=c[0], ls=ls[0], label="Radio")
ax.plot(td, FnuOGA, color=c[0], ls=ls[1], label="Optical")
ax.plot(td, FnuXGA, color=c[0], ls=ls[2], label="X-Ray")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
ax.legend()
fig.savefig("jets_gaussian.png")
plt.close(fig)

print("Plot GaussianCOre")
fig, ax = plt.subplots(1,1)
ax.plot(td, FnuRGC, color=c[0], ls=ls[0], label="Radio")
ax.plot(td, FnuOGC, color=c[0], ls=ls[1], label="Optical")
ax.plot(td, FnuXGC, color=c[0], ls=ls[2], label="X-Ray")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
ax.legend()
fig.savefig("jets_gaussiancore.png")
plt.close(fig)

print("Plot Powerlaw")
fig, ax = plt.subplots(1,1)
ax.plot(td, FnuRPL, color=c[0], ls=ls[0], label="Radio")
ax.plot(td, FnuOPL, color=c[0], ls=ls[1], label="Optical")
ax.plot(td, FnuXPL, color=c[0], ls=ls[2], label="X-Ray")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
ax.legend()
fig.savefig("jets_powerlaw.png")
plt.close(fig)

print("Plot PowerlawSmooth")
fig, ax = plt.subplots(1,1)
ax.plot(td, FnuRPS, color=c[0], ls=ls[0], label="Radio")
ax.plot(td, FnuOPS, color=c[0], ls=ls[1], label="Optical")
ax.plot(td, FnuXPS, color=c[0], ls=ls[2], label="X-Ray")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
ax.legend()
fig.savefig("jets_powerlawsmooth.png")
plt.close(fig)
