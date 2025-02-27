import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

# Set up a basic afterglow
Z = dict(jetType=grb.jet.TopHat,
         specType=grb.jet.SimpleSpec,
         thetaObs=0.0,
         E0=1.0e53,
         thetaCore=0.1,
         n0=1.0,
         p=2.5,
         epsilon_e=0.1,
         epsilon_B=0.01,
         xi_N=1.0,
         d_L=2.0e28,
         z=1.0)

# Make a second argument dict, with epsilon_e bar activated
# The "epsilon_e" parameter is now treated as "epsilon_e_bar"
#
# Join specTypes with " | "
Z_ebar = Z.copy()
Z_ebar['specType'] = grb.jet.SimpleSpec | grb.jet.EpsEBar

# Make an SED at 1 day
t = 1.0 * grb.day2sec
nu = np.geomspace(1.0e6, 1.0e20, 200)

# Compute the flux.
Fnu_p25 = grb.fluxDensity(t, nu, **Z)             # epsilon_e = 0.1
Fnu_ebar1_p25 = grb.fluxDensity(t, nu, **Z_ebar)  # epsilon_e_bar = 0.1

# Compute the epsilon_e_bar value that corresponds to the actual epsilon_e
# at p=2.5
eps_e = Z['epsilon_e']
p = Z['p']
eps_e_bar = eps_e * (p-2) / (p-1)

# Update the args and compute the flux
Z_ebar['epsilon_e'] = eps_e_bar
Fnu_ebar_p25 = grb.fluxDensity(t, nu, **Z_ebar)

# Now compute both models with p = 2.01  (not updating epsilon_e
# or epsilon_e_bar)
Z['p'] = 2.01
Z_ebar['p'] = 2.01
Fnu_p201 = grb.fluxDensity(t, nu, **Z)
Fnu_ebar_p201 = grb.fluxDensity(t, nu, **Z_ebar)

# Now just for fun, compute flux for p < 2.0.  The standard model will
# return an error, so only "ebar" here

Z_ebar['p'] = 1.5
Fnu_ebar_p15 = grb.fluxDensity(t, nu, **Z_ebar)

fig, ax = plt.subplots(1, 1)

ax.plot(nu, Fnu_p25, ls='--', color='C0', lw=2,
        label=r'Vanilla $\epsilon_e=0.1$ $p=2.5$')
ax.plot(nu, Fnu_ebar1_p25, ls='-.', color='C0', lw=1,
        label=r'Ebar $\bar{\epsilon}_e=0.1$ $p=2.5$')
ax.plot(nu, Fnu_ebar_p25, ls='-', color='C0', lw=1,
        label=r'Ebar $\epsilon_e=0.1$ $p=2.5$')

ax.plot(nu, Fnu_p201, ls='--', color='C1', lw=2,
        label=r'Vanilla $\epsilon_e=0.1$ $p=2.01$')
ax.plot(nu, Fnu_ebar_p201, ls='-', color='C1', lw=1,
        label=r'Ebar $\epsilon_e=0.1$ $p=2.01$')

ax.plot(nu, Fnu_ebar_p15, ls='-', color='C2', lw=1,
        label=r'Ebar $\epsilon_e=0.1$ $p=1.5$')

ax.set(xlabel=r'$\nu$ (Hz)', xscale='log',
       ylabel=r'$F_\nu$ (mJy)', yscale='log')

ax.legend()

figname = "example_ebar.png"
print("Saving", figname)
fig.savefig(figname)

plt.show()
