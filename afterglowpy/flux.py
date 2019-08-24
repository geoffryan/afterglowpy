from . import cocoon
from . import jet
# import time
import numpy as np


def fluxDensity(t, nu, jetType, specType, *args, **kwargs):
    """
    Compute the flux density F_nu of a GRB afterglow.

    Utiliizes the single shell approximation described in Ryan et al 2019
    to compute the synchrotron emission from the forward shock of a blast
    wave at the specified times and frequencies.

    Parameters
    ----------
    t: array_like
        Time since burst in observer frame, measured in seconds.
    nu: array_like
        Frequency of flux in observer frame, measured in Hz, same size as t.
    jetType: int
        Code for type of jet. Options are: -2 cone, -1 top hat, 0 Gaussian,
        3 quasi-spherical refreshed shock, 4 Power law.
    specType: int
        Code for type of spectrum.  Options are: 0 broken power law
        (Ryan+ 2019), 1 broken power law + inverse Compton.
    thetaObs: float
        Viewing angle in radians. umax, maximum 4-velocity of outflow, if
        jetType = 3.
    E0: float
        Isotropic-equivalent energy along the jet axis in ergs. umin, minimum
        4-velocity of outflow, if jetType = 3.
    thetaCore: float
        Half opening angle of jet core in radians. Ei, normalization of
        outflow's energy distribution in ergs, if jetType = 3.
    thetaWing: float
        Truncation angle of the jet in radians. k, power law index of
        outflow's energy distribution, if jetType = 3.
    b: float
        Power law index of jet angular energy distribution, only used if
        jetType = 4. Mej_solar, mass of material at u_max in solar masses,
        if jetType = 3.
    L0: float
        Luminosity of energy injection, in erg/s.
    q: float
        Power law index of energy injection: L = L0 (t/t0)^{-q}, t0 = 1 ks.
    ts: float
        Time energy injection ends in burster frame, in seconds.
    n0: float
        Number density of protons in circumburst medium in cm^{-3}.
    p: float
        Power law index of relativistic electron energy distribution, p > 2.
    epsilon_e: float
        Fraction of thermal energy in relativistic electrons, epsilon_e < 1.
    epsilon_B: float
        Fraction of thermal energy in magnetic field, epsilon_B < 1.
    ksiN: float
        Fraction of electrons that get accelerated, ksiN < 1.
    dL: float
        Luminosity distance to burst, in cm.
    g0: float, optional
        Initial Lorentz factor of outflow along jet axis, defaults to inf.
    z: redshift, optional
        Redshift of burst, defaults to 0.

    Other Parameters
    ----------------

    tRes: int, optional
        Time resolution, number of points per decade in t, for shock evolution.
        Defaults to 1000.
    latRes: int, optional
        Lateral resolution of structured jet, number of conical sections per
        thetaCore-sized interval. Defaults to 5.
    rtol: float, optional
        Relative tolerance of flux integration, defaults to 1.0e-4.
    spread: {'True', 'False'}
        Whether to include jet spreading. Defaults to True.

    Returns
    -------

    The flux density F_nu in the observer frame, same shape as t and nu.
    """

    if 'z' in kwargs.keys():
        z = kwargs.pop('z')
    else:
        z = 0.0

    tz = t / (1+z)
    nuz = nu * (1+z)

    # Default spreading method
    if 'spread' in kwargs:
        if kwargs['spread'] is True:
            if jetType == -2 and 'thetaCoreGlobal' in kwargs:
                kwargs['spread'] = 8
            else:
                kwargs['spread'] = 7

    # timeA = time.time()
    if jetType == 3:
        Fnu = cocoon.fluxDensity(tz, nuz, jetType, specType, *args, **kwargs)
    else:
        Fnu = jet.fluxDensity(tz, nuz, jetType, specType, *args, **kwargs)
    # timeB = time.time()
    # print("Eval took: {0:f} s".format(timeB - timeA))

    # K-correct the flux
    Fnu *= 1+z

    return Fnu


def intensity(theta, phi, t, nu, jetType, specType, *args, **kwargs):
    """
    Compute the specific intensity I_nu of a GRB afterglow.

    Utiliizes the single shell approximation described in Ryan et al 2019
    to compute the synchrotron emission from the forward shock of a blast
    wave at the specified angular coordinates, times, and frequencies.

    Angular coordinates are in a spherical coordinate system, centered on the
    burst, with z-axis aligned on the jet axis.

    Parameters
    ----------
    theta: array_like
        Polar angle from jet axis in radians.
    phi: array_like
        Azimuthal angle around jet axis in radians. Observer is at phi = 0.
        Same size as theta.
    t: array_like
        Time since burst in observer frame, measured in seconds. Same size as
        theta.
    nu: array_like
        Frequency of flux in observer frame, measured in Hz, same size as
        theta.
    jetType: int
        Code for type of jet. Options are: -2 cone, -1 top hat, 0 Gaussian,
        3 quasi-spherical refreshed shock, 4 Power law.
    specType: int
        Code for type of spectrum.  Options are: 0 broken power law
        (Ryan+ 2019), 1 broken power law + inverse Compton.
    thetaObs: float
        Viewing angle in radians. umax, maximum 4-velocity of outflow, if
        jetType = 3.
    E0: float
        Isotropic-equivalent energy along the jet axis in ergs. umin, minimum
        4-velocity of outflow, if jetType = 3.
    thetaCore: float
        Half opening angle of jet core in radians. Ei, normalization of
        outflow's energy distribution in ergs, if jetType = 3.
    thetaWing: float
        Truncation angle of the jet in radians. k, power law index of
        outflow's energy distribution, if jetType = 3.
    b: float
        Power law index of jet angular energy distribution, only used if
        jetType = 4. Mej_solar, mass of material at u_max in solar masses,
        if jetType = 3.
    L0: float
        Luminosity of energy injection, in erg/s.
    q: float
        Power law index of energy injection: L = L0 (t/t0)^{-q}, t0 = 1 ks.
    ts: float
        Time energy injection ends in burster frame, in seconds.
    n0: float
        Number density of protons in circumburst medium in cm^{-3}.
    p: float
        Power law index of relativistic electron energy distribution, p > 2.
    epsilon_e: float
        Fraction of thermal energy in relativistic electrons, epsilon_e < 1.
    epsilon_B: float
        Fraction of thermal energy in magnetic field, epsilon_B < 1.
    ksiN: float
        Fraction of electrons that get accelerated, ksiN < 1.
    dL: float
        Luminosity distance to burst, in cm.
    g0: float, optional
        Initial Lorentz factor of outflow along jet axis, defaults to inf.
    z: redshift, optional
        Redshift of burst, defaults to 0.

    Other Parameters
    ----------------

    tRes: int, optional
        Time resolution, number of points per decade in t, for shock evolution.
        Defaults to 1000.
    latRes: int, optional
        Lateral resolution of structured jet, number of conical sections per
        thetaCore-sized interval. Defaults to 5.
    rtol: float, optional
        Relative tolerance of flux integration, defaults to 1.0e-4.
    spread: {'True', 'False'}
        Whether to include jet spreading. Defaults to True.

    Returns
    -------

    The flux density F_nu in the observer frame, same shape as t and nu.
    """

    if 'z' in kwargs.keys():
        z = kwargs.pop('z')
    else:
        z = 0.0

    tz = t / (1+z)
    nuz = nu * (1+z)

    # Default spreading method
    if 'spread' in kwargs:
        if kwargs['spread'] is True:
            if jetType == -2 and 'thetaCoreGlobal' in kwargs:
                kwargs['spread'] = 8
            else:
                kwargs['spread'] = 7

    Inu = np.empty(theta.shape)
    Inu.flat[:] = jet.intensity(theta.flat, phi.flat, tz.flat, nuz.flat,
                                jetType, specType, *args, **kwargs)

    # K-correct the intensity
    # I'm only using the flux correction here, which leaves the angular
    # part of the intensity uncorrected.  Best be careful.

    Inu *= 1+z

    return Inu
