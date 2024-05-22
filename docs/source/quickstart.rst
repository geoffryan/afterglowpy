Quickstart
==========

*afterglowpy* exposes a single function ``fluxDensity()`` for the purposes of calculating light curves and spectra from GRB afterglows.  It takes many arguments, which are grouped as follows:

1. ``t``: a 1-D array of observer-frame times, in seconds (s)
2. ``nu``: a 1-D array, the same shape as ``t``, of observer-frame frequencies, in Hertz (Hz)
3. ``jetType``: an integer specifying which jet model to use (top-hat: -1, Gaussian: 0, power law: 4).
4. ``specType``: an integer specifying which spectral model to use (must be 0).
5. Model Parameters (keywords): A model-dependent set of parameters, such as E\ :sub:`iso`, n\ :sub:`0`, etc.
6. Numerical Parameters (keywords): A model-dependent set of keyword arguments specifying numerical settings, such as `tRes` (resolution of time integration), and `spread` (whether or not to include jet spreading).

A light curve from a top-hat jet may be created with::
    
    import numpy as np
    import afterglowpy as grb

    # For convenience, place arguments into a dict.
    Z = {'jetType':     grb.jet.TopHat,     # Top-Hat jet
         'specType':    grb.jet.SimpleSpec, # Basic Synchrotron Spectrum

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

``Fnu`` is now an array, same size as ``t`` and ``nu``, containing the observed flux in mJy at each time and frequency.

Let's calculate the light curve for a Gaussian jet with the same parameters. A Gaussian jet has an isotropic-energy profile E(theta) = E0 * exp(-0.5*theta^2/thetaCore^2).  Now that ``thetaCore`` sets the Gaussian width of the jet core, and we need to provide a sensible value for ``thetaWing`` to truncate the outer regions of the jet::

    Z['jetType'] = grb.jet.Gaussian     # Gaussian jet

    # We'll re-use the Z dict, just add a thetaWing entry
    
    Z['thetaObs'] = 0.3                  # Larger viewing angle, just for fun
    Z['thetaWing'] = 4 * Z['thetaCore']  # Setting thetaW

    # Nothing else to change, so calculate!

    Fnu_G = grb.fluxDensity(t, nu, **Z)

A power law structured jet has E(theta) ~ theta^(-b) for theta >> thetaC.  We need to give ``b`` a reasonable value::

    Z['jetType'] = grb.jet.PowerLaw     # power law

    # Set b
    Z['b'] = 6.0      # The rough ballpark from RHD simulations

    # Calculate!

    Fnu_PL = grb.fluxDensity(t, nu, **Z)

There you go! Simple X-ray light curves for three different jet models.
