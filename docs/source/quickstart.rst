Quickstart
==========

*afterglowpy* exposes a single function ``fluxDensity()`` for the purposes of calculating light curves and spectra from GRB afterglows.  It takes many arguments, which are grouped as follows:

1. ``t``: a 1-D array of observer-frame times, in seconds (s)
2. ``nu``: a 1-D array, the same shape as ``t``, of observer-frame frequencies, in Hertz (Hz)
3. ``jetType``: an integer specifying which jet model to use (top-hat: -1, Gaussian: 0, power law: 4).
4. ``specType``: an integer specifying which spectral model to use (must be 0).
5. Model Parameters (positional): A model-dependent set of parameters, such as E\ :sub:`iso`, n\ :sub:`0`, etc.
6. Numerical Parameters (keywords): A model-dependent set of keyword arguments specifying numerical settings, such as `tRes` (resolution of time integration), and the redshift.

A light curve from a top-hat jet may be created with::
    
    import numpy as np
    import afterglowpy as grb

    jetType = -1     # top hat jet

    thetaObs = 0.05  # Viewing angle in radians
    E0 = 1.0e53      # Isotropic-equivalent energy in erg
    thetaC = 0.1     # Half-opening angle in radians
    thetaW = 0.1     # Truncation angle, unused for top-hat
    b = 0            # power law index, unused for top-hat
    n0 = 1.0         # circumburst density in cm^{-3}
    p = 2.2          # electron energy distribution index
    eps_e = 0.1      # epsilon_e
    eps_B = 0.01     # epsilon_B
    xi_N = 1.0       # Fraction of electrons accelerated
    dL = 1.0e28      # Luminosity distance in cm
    z = 0.55         # redshift

    # Space time points geometrically, from 10^3 s to 10^7 s
    t = np.geomspace(1.0e3, 1.0e7, 300)

    # Calculate flux in a single X-ray band (all times have same frequency)
    nu = np.empty(t.shape)
    nu[:] = 1.0e18

    # For convenience, place positional arguments in an array, 
    # and keywords into a dict

    Y = np.array([thetaObs, E0, thetaC, thetaW, b, 0, 0, 0,
                  n0, p, epse_e, eps_B, xi_N, dL])
    Z = {'z': z}

    # Calculate!

    Fnu = grb.fluxDensity(t, nu, jetType, 0, *Y, **Z)

``Fnu`` is now an array, same size as ``t`` and ``nu``, containing the observed flux in mJy at each time and frequency.

Let's calculate the light curve for a Gaussian jet with the same parameters. A Gaussian jet has an isotropic-energy profile E(theta) = E0 * exp(-0.5*theta^2/thetaC^2).  Now that ``thetaC`` sets the Gaussian width of the jet core, and we need to provide a sensible value for ``thetaW`` to truncate the outer regions of the jet::

    jetType = 0     # Gaussian jet

    # We'll re-use the Y[] array, just change some values in place
    
    Y[0] = 0.3       # Larger viewing angle, just for fun
    Y[3] = 4*thetaC  # Setting thetaW

    # Nothing else to change, so calculate!

    Fnu_G = grb.fluxDensity(t, nu, jetType, 0, *Y, **Z)

A power law structured jet has E(theta) ~ theta^(-b) for theta >> thetaC.  We need to give ``b`` a reasonable value::

    jetType = 4     # power law

    # Set b
    Y[4] = 6.0      # The rough ballpark from RHD simulations

    # Calculate!

    Fnu_PL = grb.fluxDensity(t, nu, jetType, 0, *Y, **Z)

There you go! Simple X-ray light curves for three different jet models.
