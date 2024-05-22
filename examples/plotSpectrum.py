import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

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


nua = 1.0e0   # Low Frequencies in Hz
nub = 1.0e20  # High Frequencies in Hz

t = 1.0 * grb.day2sec  # spectrum at 1 day
nu = np.geomspace(nua, nub, num=100)

print("Calculating")
Fnu = grb.fluxDensity(t, nu, **Z)

print("Writing spec.txt")
with open("spec.txt", 'w') as f:
    f.write("# t " + str(t) + ' (s)\n')
    f.write("# nu(Hz)   Fnu(mJy)\n")
    for i in range(len(nu)):
        f.write("{0:.6e} {1:.6e}\n".format(nu[i], Fnu[i]))

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
