import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

# For convenience, place arguments into a dict.
Z = {'jetType':     grb.jet.TopHat,     # Top-Hat jet
     'specType':    grb.jet.SimpleSpec, # Basic Synchrotron Emission Spectrum

     'thetaObs':    0.05,   # Viewing angle in radians
     'E0':          1.0e53, # Isotropic-equivalent energy in erg
     'thetaCore':   0.1,    # Half-opening angle in radians
     'n0':          1.0,    # circumburst density in cm^{-3}
     'p':           2.2,    # electron energy distribution index
     'epsilon_e':   0.1,    # epsilon_e
     'epsilon_B':   0.01,   # epsilon_B
     'xi_N':        1.0,    # Fraction of electrons accelerated
     'd_L':         1.0e28, # Luminosity distance in cm
     'z':           0.55}   # redshift

# Space time points geometrically, from 10^3 s to 10^7 s
t = np.geomspace(1.0e3, 1.0e7, 300)

# Calculate flux in a single X-ray band (all times have same frequency)
nu = np.empty(t.shape)
nu[:] = 1.0e18

# Calculate!

Fnu = grb.fluxDensity(t, nu, **Z)

# Write to a file

print("Writing lc.txt")
with open("lc.txt", 'w') as f:
    f.write("# nu " + str(nu) + '\n')
    f.write("# t(s)     Fnu(mJy)\n")
    for i in range(len(t)):
        f.write("{0:.6e} {1:.6e}\n".format(t[i], Fnu[i]))

# Plot!

print("Plotting")
fig, ax = plt.subplots(1, 1)

ax.plot(t, Fnu)

ax.set(xscale='log', xlabel=r'$t$ (s)',
       yscale='log', ylabel=r'$F_\nu$[$10^{18}$ Hz] (mJy)')

fig.tight_layout()
print("Saving figure lc.png")
fig.savefig("lc.png")
plt.close(fig)
