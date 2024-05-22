import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

def computeMoments(t, nu, Z):

    # Make array to store moment options
    moment = np.empty(t.shape, dtype=int)

    # compute moment 0 (flux)
    moment[:] = grb.jet.MOM_0
    Fnu = grb.fluxDensity(t, nu, **Z, moment=moment)

    # compute moment 1 (integral of x_obs * intensity)
    moment[:] = grb.jet.MOM_X
    FnuX = grb.fluxDensity(t, nu, **Z, moment=moment)

    # compute moment 2 (integral of x_obs^2 * intensity)
    moment[:] = grb.jet.MOM_XX
    FnuXX = grb.fluxDensity(t, nu, **Z, moment=moment)

    # compute moment 2 (integral of y_obs^2 * intensity)
    moment[:] = grb.jet.MOM_YY
    FnuYY = grb.fluxDensity(t, nu, **Z, moment=moment)

    # The jet propagates in the x-direction and is axisymmetric, so MOM_Y
    # and MOM_XY will always be identically 0

    # Get the intensity-weighted distance measures
    X_cm = FnuX / Fnu     # in cm
    X2_cm = FnuXX / Fnu   # in cm^2
    Y2_cm = FnuYY / Fnu   # in cm^2

    # Convert proper distances to angular distances
    dA = Z['d_L'] / (1 + Z['z'])**2

    X_rad = X_cm / dA
    X2_rad = X2_cm / dA**2
    Y2_rad = Y2_cm / dA**2

    # Convert radians to mas
    rad2mas = 1000 * 3600 * 180 / np.pi
    x = X_rad * rad2mas
    x2 = X2_rad * rad2mas**2
    y2 = Y2_rad * rad2mas**2

    # Compute sizes
    sig_x = np.sqrt(x2 - x**2)
    sig_y = np.sqrt(y2)

    return Fnu, x, sig_x, sig_y

if __name__ == "__main__":


# For convenience, place arguments into a dict.
    Z = {'jetType':     grb.jet.Gaussian,     # Top-Hat jet
         'specType':    grb.jet.SimpleSpec, # Basic Synchrotron Spectrum

         'thetaObs':    0.00,   # Viewing angle in radians
         'E0':          1.0e53, # Isotropic-equivalent energy in erg
         'thetaCore':   0.05,    # Half-opening angle in radians
         'thetaWing':   0.4,    # Half-opening angle in radians
         'n0':          1.0e-3,    # circumburst density in cm^{-3}
         'p':           2.2,    # electron energy distribution index
         'epsilon_e':   0.1,    # epsilon_e
         'epsilon_B':   0.01,   # epsilon_B
         'xi_N':        1.0,    # Fraction of electrons accelerated
         'd_L':         1.0e26, # Luminosity distance in cm
         'z':           0.01}   # redshift

    # Space time points geometrically, from 10^3 s to 10^7 s
    t = np.geomspace(1.0e3, 1.0e7, 100)

    # Calculate flux in a single X-ray band (all times have same frequency)
    nu = np.empty(t.shape)
    nu[:] = 3.0e9


    # Calculate!
    Fnu, x, sig_x, sig_y = computeMoments(t, nu, Z)
    
    # Try an off-axis one too
    Z['thetaObs'] = 6 * Z['thetaCore']
    Fnu_off, x_off, sig_x_off, sig_y_off = computeMoments(t, nu, Z)

    # Plot!

    print("Plotting")
    fig, ax = plt.subplots(4, 1, figsize=(6, 12))

    ax[0].plot(t, Fnu, label='on-axis')
    ax[0].plot(t, Fnu_off, label='off-axis')
    
    ax[1].plot(t, x)
    ax[1].plot(t, x_off)
    
    ax[2].plot(t, sig_x)
    ax[2].plot(t, sig_x_off)
    
    ax[3].plot(t, sig_y)
    ax[3].plot(t, sig_y_off)

    ax[0].legend()

    ax[0].set(xscale='log', xlabel=r'$t$ (s)',
              yscale='log', ylabel=r'$F_\nu$[$3\times 10^{9}$ Hz] (mJy)')
    ax[1].set(xscale='log', xlabel=r'$t$ (s)',
              yscale='linear', ylabel=r'centroid offset (mas)')
    ax[2].set(xscale='log', xlabel=r'$t$ (s)',
              yscale='log', ylabel=r'$\sigma_x$ (mas)')
    ax[3].set(xscale='log', xlabel=r'$t$ (s)',
              yscale='log', ylabel=r'$\sigma_y$ (mas)')

    fig.tight_layout()
    print("Saving figure lc_centroid_and_size.png")
    fig.savefig("lc_centroid_and_size.png")
    plt.close(fig)
