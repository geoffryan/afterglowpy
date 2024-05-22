import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

# Jet Parameters
Z = {'jetType':     grb.jet.Gaussian,     # Gaussian jet
     'specType':    grb.jet.SimpleSpec,   # Basic Synchrotron Emission Spectrum

     'thetaObs':    0.3,   # Viewing angle in radians
     'E0':          1.0e53, # Isotropic-equivalent energy in erg
     'thetaCore':   0.05,    # Half-opening angle in radians
     'thetaWing':   0.4,    # Outer truncation angle
     'n0':          1.0e-3,    # circumburst density in cm^{-3}
     'p':           2.2,    # electron energy distribution index
     'epsilon_e':   0.1,    # epsilon_e
     'epsilon_B':   0.0001,   # epsilon_B
     'xi_N':        1.0,    # Fraction of electrons accelerated
     'd_L':         1.36e26, # Luminosity distance in cm
     'z':           0.01}   # redshift

# Time and Frequencies
ta = 1.0e-1 * grb.day2sec
tb = 1.0e3 * grb.day2sec
t = np.geomspace(ta, tb, num=100)

nuR = 6.0e9
nuO = 1.0e14
nuX = 1.0e18

# Calculate!
print("Calc Radio")
FnuR = grb.fluxDensity(t, nuR, **Z)
print("Calc Optical")
FnuO = grb.fluxDensity(t, nuO, **Z)
print("Calc X-ray")
FnuX = grb.fluxDensity(t, nuX, **Z)

# Plot!
print("Plot")

tday = t * grb.sec2day

fig, ax = plt.subplots(1, 1)
ax.plot(tday, FnuR, ls='-', label=r'$\nu=6$ GHz')
ax.plot(tday, FnuO, ls='--', label=r'$\nu=10^{14}$ Hz')
ax.plot(tday, FnuX, ls='-.', label=r'$\nu=10^{18}$ Hz')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$t$ (d)')
ax.set_ylabel(r'$F_\nu$ (mJy)')
ax.legend()
fig.tight_layout()

print("Saving lc_multi.png")
fig.savefig("lc_multi.png")
plt.close(fig)
