import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

# For convenience, place arguments into a dict.
Z = {'jetType':     grb.jet.TopHat,     # Top-Hat jet
     'specType':    grb.jet.DeepNewtonian, # Deep Newtonian Emission

     'thetaObs':    0.05,   # Viewing angle in radians
     'E0':          1.0e53, # Isotropic-equivalent energy in erg
     'thetaCore':   0.1,    # Half-opening angle in radians
     'n0':          1.0,    # circumburst density in cm^{-3}
     'p':           2.2,    # electron energy distribution index
     'epsilon_e':   0.001,  # epsilon_e (low eps_e turns DN on sooner)
     'epsilon_B':   0.001,  # epsilon_B
     'xi_N':        1.0,    # Fraction of electrons accelerated
     'd_L':         1.0e28, # Luminosity distance in cm
     'z':           0.55,   # redshift
     'counterjet':  True    # Turn on the counter jet
     }

# Space time points geometrically, from 10^4 s to 10^9 s
t = np.geomspace(1.0e4, 1.0e9, 300)

# Calculate flux in a single radio band (all times have same frequency)
nu = np.empty(t.shape)
nu[:] = 3.0e9

# Calculate!

Fnu = grb.fluxDensity(t, nu, **Z)

# Compute a light curve with counterjet but no deep newtonian

Z_cj_only = Z.copy()
Z_cj_only['specType'] = grb.jet.SimpleSpec
Fnu_cj_only = grb.fluxDensity(t, nu, **Z_cj_only)

# Compute a light curve with no counterjet and no deep newtonian
Z_default = Z.copy()
Z_default['specType'] = grb.jet.SimpleSpec
Z_default['counterjet'] = False
Fnu_default = grb.fluxDensity(t, nu, **Z_default)


# Plot!

print("Plotting")
fig, ax = plt.subplots(1, 1)

ax.plot(t, Fnu, label='Deep Newtonian and counterjet')
ax.plot(t, Fnu_cj_only, label='counterjet')
ax.plot(t, Fnu_default, label='Default (no Deep Newtonian or counterjet)')
ax.legend()

ax.set(xscale='log', xlabel=r'$t$ (s)',
       yscale='log', ylabel=r'$F_\nu$[$3\times 10^{9}$ Hz] (mJy)')

fig.tight_layout()
print("Saving figure lc_late_time.png")
fig.savefig("lc_late_time.png")
plt.close(fig)
