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
    t: array_like or scalar
        Time since burst in observer frame, measured in seconds.
    nu: array_like or scalar
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

    Raises
    ------

    ValueError
        If t, nu are the wrong shape or arguments take illegal values.
    """

    # Check Arguments, will raise ValueError if args are bad
    t, nu = checkTNu(t, nu)

    if jetType == 3:
        checkCocoonArgs(jetType, specType, *args, **kwargs)
    else:
        checkJetArgs(jetType, specType, *args, **kwargs)

    # arguments are good, full steam ahead!
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

    Fnu = np.empty(tz.shape)

    if jetType == 3:
        Fnu.flat[:] = cocoon.fluxDensity(tz.flat, nuz.flat, jetType, specType,
                                         *args, **kwargs)
    else:
        Fnu.flat[:] = jet.fluxDensity(tz.flat, nuz.flat, jetType, specType,
                                      *args, **kwargs)
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
    theta: array_like or scalar
        Polar angle from jet axis in radians. Scalar, or array of same shape
        as phi, t, and nu.
    phi: array_like or scalar
        Azimuthal angle around jet axis in radians. Observer is at phi = 0.
        Scalar, or array of same shape as theta, t, and nu.
    t: array_like or scalar
        Time since burst in observer frame, measured in seconds. Scalar, or
        array of same shape as theta, phi, and nu.
    nu: array_like or scalar
        Frequency of flux in observer frame, measured in Hz. Scalar, or array
        of same shape as theta, phi, and t.
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

    Raises
    ------

    ValueError
        If theta, phi, t, nu are the wrong shape or arguments take illegal
        values.
    """

    # Check Arguments, will raise ValueError if args are bad
    t, nu = checkThetaPhiTNu(theta, phi, t, nu)

    if jetType == 3:
        checkCocoonArgs(jetType, specType, *args, **kwargs)
    else:
        checkJetArgs(jetType, specType, *args, **kwargs)

    # arguments are good, full steam ahead!

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


def checkTNu(t, nu):
    # Make sure t and nu are array_like or castable to an array.
    t = np.atleast_1d(t)
    nu = np.atleast_1d(nu)

    # Check shapes, if scalars make into right size array
    if t.shape != nu.shape:
        if t.shape == (1, ):
            T = t[0]
            t = np.empty(nu.shape)
            t[:] = T
        elif nu.shape == (1, ):
            NU = nu[0]
            nu = np.empty(t.shape)
            nu[:] = NU
        else:
            raise ValueError("t and nu must be same shape or scalars")

    return t, nu


def checkThetaPhiTNu(theta, phi, t, nu):
    # Make sure args are array_like or castable to an array.
    theta = np.atleast_1d(theta)
    phi = np.atleast_1d(phi)
    t = np.atleast_1d(t)
    nu = np.atleast_1d(nu)

    shape = (1, )
    if shape != theta.shape:
        shape = theta.shape
    elif shape != phi.shape:
        shape = phi.shape
    elif shape != t.shape:
        shape = t.shape
    else:
        shape = nu.shape

    if theta.shape != shape:
        if theta.shape == (1, ):
            TH = theta[0]
            theta = np.empty(shape)
            theta[:] = TH
        else:
            msg = "theta must be scalar or same shape as phi, t, nu"
            raise ValueError(msg)

    if phi.shape != shape:
        if phi.shape == (1, ):
            PH = phi[0]
            phi = np.empty(shape)
            phi[:] = PH
        else:
            msg = "phi must be scalar or same shape as theta, t, nu"
            raise ValueError(msg)

    if t.shape != shape:
        if t.shape == (1, ):
            T = t[0]
            t = np.empty(shape)
            t[:] = T
        else:
            msg = "t must be scalar or same shape as theta, phi, nu"
            raise ValueError(msg)

    if nu.shape != shape:
        if nu.shape == (1, ):
            NU = nu[0]
            nu = np.empty(shape)
            nu[:] = NU
        else:
            msg = "nu must be scalar or same shape as theta, phi, t"
            raise ValueError(msg)

    return theta, phi, t, nu


def checkJetArgs(jetType, specType, *args, **kwargs):

    for x in args:
        if not np.isfinite(x):
            raise ValueError("All parameters must be finite")

    theta_obs = args[0]
    E0 = args[1]
    theta_c = args[2]
    theta_w = args[3]
    b = args[4]
    L0 = args[5]
    t_s = args[7]
    n0 = args[8]
    p = args[9]
    epse = args[10]
    epsB = args[11]
    xiN = args[12]
    dL = args[13]

    # More-or-less universal bounds
    if theta_obs < 0.0 or theta_obs > 0.5*np.pi:
        raise ValueError("theta_obs must be in [0.0, pi/2]")
    if E0 <= 0.0:
        raise ValueError("E0 must be positive")
    if theta_c <= 0.0 or theta_c > 0.5*np.pi:
        raise ValueError("theta_c must be in (0.0, pi/2]")
    if jetType != -1 and (theta_w <= 0.0 or theta_w > 0.5*np.pi):
        raise ValueError("theta_w must be in (0.0, pi/2]")
    if jetType == 4 and b <= 0.0:
        raise ValueError("b must be positive")
    if L0 < 0.0:
        raise ValueError("L0 must be non-negative")
    if t_s < 0.0:
        raise ValueError("t_s must be non-negative")
    if n0 <= 0.0:
        raise ValueError("n0 must be positive")
    if specType != 2 and p <= 2.0:
        raise ValueError("p must be in (2, inf)")
    if specType == 2 and p <= 1.0:
        raise ValueError("p must be in (1, inf)")
    if epse <= 0.0 or epse > 1.0:
        raise ValueError("epsilon_e must be in (0, 1]")
    if epsB <= 0.0 or epsB > 1.0:
        raise ValueError("epsilon_B must be in (0, 1]")
    if xiN <= 0.0 or xiN > 1.0:
        raise ValueError("xi_N must be in (0, 1]")
    if dL <= 0.0:
        raise ValueError("dL must be positive")

    if 'z' in kwargs:
        if kwargs['z'] < 0.0:
            raise ValueError("z must be non-negative")

    # Model Specific bounds
    if jetType == -2 and theta_c > theta_w:
        raise ValueError("theta_w must be larger than theta_c for cone model")

    return


def checkCocoonArgs(jetType, specType, *args, **kwargs):

    for x in args:
        if not np.isfinite(x):
            raise ValueError("All parameters must be finite")

    u_max = args[0]
    u_min = args[1]
    Ei = args[2]
    Mej_solar = args[4]
    L0 = args[5]
    t_s = args[7]
    n0 = args[8]
    p = args[9]
    epse = args[10]
    epsB = args[11]
    xiN = args[12]
    dL = args[13]

    if u_max <= 0.0:
        raise ValueError("u_max must be positive")
    if u_min <= 0.0:
        raise ValueError("u_min must be positive")
    if Ei <= 0.0:
        raise ValueError("Ei must be positive")
    if Mej_solar <= 0.0:
        raise ValueError("Mej_solar must be positive")
    if L0 < 0.0:
        raise ValueError("L0 must be non-negative")
    if t_s < 0.0:
        raise ValueError("t_s must be non-negative")
    if n0 <= 0.0:
        raise ValueError("n0 must be positive")
    if specType != 2 and p <= 2.0:
        raise ValueError("p must be in (2, inf)")
    if epse <= 0.0 or epse > 1.0:
        raise ValueError("epsilon_e must be in (0, 1]")
    if epsB <= 0.0 or epsB > 1.0:
        raise ValueError("epsilon_B must be in (0, 1]")
    if xiN <= 0.0 or xiN > 1.0:
        raise ValueError("xi_N must be in (0, 1]")
    if dL <= 0.0:
        raise ValueError("dL must be positive")

    if 'z' in kwargs:
        if kwargs['z'] < 0.0:
            raise ValueError("z must be non-negative")

    return
