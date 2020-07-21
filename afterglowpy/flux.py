from . import cocoon
from . import jet
import numpy as np


JET_DEFAULT = {
    "thetaObs":         0.0,
    "E0":               1.0e53,
    "thetaCore":        0.1,
    "thetaWing":        0.4,
    "b":                4.0,
    "L0":               0.0,
    "q":                0.0,
    "ts":               0.0,
    "n0":               1.0,
    "p":                2.2,
    "epsilon_e":        0.1,
    "epsilon_B":        0.01,
    "ksiN":             1.0,
    "dL":               1.0e27,
    "g0":               1.0e3,
    "E0Global":         1.0e53,
    "thetaCoreGlobal":  0.1,
    "tRes":             1000,
    "latRes":           5,
    "rtol":             1.0e-4,
    "mask":             None,
    "spread":           7,
    "gammaType":        0
    }

JET_KEYS = list(JET_DEFAULT.keys())


def fluxDensity(t, nu, jetType=-1, specType=0,
                thetaObs=0.0, E0=1.0e53, thetaCore=0.1, thetaWing=0.4,
                b=4.0, L0=0.0, q=0.0, ts=0.0, n0=1.0, p=2.2, epsilon_e=0.1,
                epsilon_B=0.01, ksiN=1.0, dL=1.0e27, g0=-1.0, LR=0.0, LO=0.0,
                LX=0.0, tAdd=0.0, z=0.0,
                **kwargs):
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
        3 quasi-spherical refreshed shock, 4 Power law. Default: -1
    specType: int
        Code for type of spectrum.  Options are: 0 broken power law
        (Ryan+ 2019), 1 broken power law w/ inverse Compton cooling. Default: 0
    thetaObs: float
        Viewing angle in radians. umax, maximum 4-velocity of outflow, if
        jetType = 3. Default: 0.0.
    E0: float
        Isotropic-equivalent energy along the jet axis in ergs. umin, minimum
        4-velocity of outflow, if jetType = 3. Default: 1.0e53.
    thetaCore: float
        Half opening angle of jet core in radians. Ei, normalization of
        outflow's energy distribution in ergs, if jetType = 3. Default 0.1.
    thetaWing: float
        Truncation angle of the jet in radians. k, power law index of
        outflow's energy distribution, if jetType = 3. Default 0.4.
    b: float
        Power law index of jet angular energy distribution, only used if
        jetType = 4. Mej_solar, mass of material at u_max in solar masses,
        if jetType = 3.  Default: 4.
    L0: float
        Luminosity of energy injection, in erg/s.  Default 0.0.
    q: float
        Power law index of energy injection: L = L0 (t/t0)^{-q}, t0 = 1 ks.
        Default 0.0.
    ts: float
        Time energy injection ends in burster frame, in seconds. Default 0.0.
    n0: float
        Number density of protons in circumburst medium in cm^{-3}.
        Default 1.0.
    p: float
        Power law index of relativistic electron energy distribution, p > 2.
        Default 2.2.
    epsilon_e: float
        Fraction of thermal energy in relativistic electrons, epsilon_e < 1.
        Default 0.1.
    epsilon_B: float
        Fraction of thermal energy in magnetic field, epsilon_B < 1.
        Default 0.01.
    ksiN: float
        Fraction of electrons that get accelerated, ksiN < 1. Default 1.0.
    dL: float
        Luminosity distance to burst, in cm. Default 1.0e27.
    g0: float, optional
        Initial Lorentz factor of outflow along jet axis, defaults to inf.
    LR: float, optional
        Additional constant source-frame radio luminosity to add to the jet
        component.  Applied uniformly for nu < 300 GHz. L_nu
        normalized with 10 GHz bandwidth: L_nu = L / 10 GHz.
        Non-negative, erg/s.
    LO: float, optional
        Additional constant source-frame optical luminosity to add to the jet
        component.  Applied uniformly for 300 GHz < nu < 0.1 keV (~2.4e16 Hz).
        L_nu normalized with 2324.78 A bandwidth (FWHM of F606W):
        L_nu = L * 2324.78 A / c. Non-negative, erg/s.
    LX: float, optional
        Additional constant source-frame X-ray luminosity to add to the jet
        component.  Applied uniformly for nu > 0.1 keV (~2.4e16 Hz).
        L_nu normalized with 0.3 keV - 10 keV bandwidth: L_nu = L / 9.7 keV.
        Non-negative, erg/s.
    tAdd: float, optional
        Time at which additional luminosities begin, in seconds.
    z: float, optional
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
    Fnu: array
        The flux density F_nu in the observer frame, same shape as t and nu.

    Raises
    ------

    ValueError
        If t, nu are the wrong shape or arguments take illegal values.
    """

    # Package all the arguments together.
    # This would be easier if we could just use **kwargs alone, but 
    # we need to keep the ordering of the arguments for people who use
    # the positional interface

    argsDict = kwargs.copy()
    for key, val in locals().items():
        if (key != 'kwargs' and key != 't' and key != 'nu'
                and key != 'argsDict'):
            argsDict[key] = val
    
    # Check Arguments, will raise ValueError if args are bad
    t, nu = checkTNu(t, nu)

    if jetType == jet.Spherical:
        checkCocoonArgs(**argsDict)
    else:
        checkJetArgs(**argsDict)

    # arguments are good, full steam ahead!
    if 'z' in argsDict.keys():
        z = argsDict.pop('z')
    else:
        z = 0.0

    tz = t / (1+z)
    nuz = nu * (1+z)

    # Default spreading method
    if 'spread' in argsDict:
        if argsDict['spread'] == True:
            if jetType == -2 and 'thetaCoreGlobal' in argsDict:
                argsDict['spread'] = 8
            else:
                argsDict['spread'] = 7

    # Intercept background luminosities

    if 'LR' in argsDict.keys():
        argsDict.pop('LR')
    if 'LO' in argsDict.keys():
        argsDict.pop('LO')
    if 'LX' in argsDict.keys():
        argsDict.pop('LX')
    if 'tAdd' in argsDict.keys():
        argsDict.pop('tAdd')

    # timeA = time.time()

    Fnu = np.empty(tz.shape)

    if jetType == jet.Spherical:
        Fnu.flat[:] = cocoon.fluxDensity(tz.flat, nuz.flat, **argsDict)
    else:
        Fnu.flat[:] = jet.fluxDensity(tz.flat, nuz.flat, **argsDict)
    # timeB = time.time()
    # print("Eval took: {0:f} s".format(timeB - timeA))

    
    # Adding background luminosities.
    L_to_flux = cocoon.cgs2mJy / (4*np.pi*dL*dL)


    if LR > 0.0:
        rad = (nuz < 3.0e11) & (tz > tAdd)  # radio < 300 GHz
        Lnu = LR / 1.0e10  # 10 GHz bandwidth
        Fnu[rad] += Lnu*L_to_flux
    if LO > 0.0:
        # 300GHz<opt<100eV
        opt = (nuz >= 3.0e11) & (nuz < 100*cocoon.eV2Hz) & (tz > tAdd)
        Lnu = LO * 2.32478e-5 / cocoon.c # 2324.78 A bandwidth
        Fnu[opt] += Lnu*L_to_flux
    if LX > 0.0:
        xry = (nuz >= 100*cocoon.eV2Hz) & (tz > tAdd)  # xray > 100eV
        Lnu = LX / (9.7e3 * cocoon.eV2Hz)  # 9.7 keV bandwidth
        Fnu[xry] += Lnu*L_to_flux

    # K-correct the flux
    Fnu *= 1+z

    return Fnu


def intensity(theta, phi, t, nu, jetType=-1, specType=0,
              thetaObs=0.0, E0=1.0e53, thetaCore=0.1, thetaWing=0.4,
              b=4.0, L0=0.0, q=0.0, ts=0.0, n0=1.0, p=2.2, epsilon_e=0.1,
              epsilon_B=0.01, ksiN=1.0, dL=1.0e27, g0=-1.0, LR=0.0, LO=0.0,
              LX=0.0, tAdd=0.0, z=0.0,
              **kwargs):
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
    Inu : array
        The specific intensity I_nu in the observer frame, same shape as
        theta, phi, t, and nu.

    Raises
    ------

    ValueError
        If theta, phi, t, nu are the wrong shape or arguments take illegal
        values.
    """
    
    argsDict = kwargs.copy()
    for key, val in locals().items():
        if (key != 'kwargs' and key != 't' and key != 'nu'
                and key != 'theta' and key != 'phi'
                and key != 'argsDict'):
            argsDict[key] = val

    # Check Arguments, will raise ValueError if args are bad
    theta, phi, t, nu = checkThetaPhiTNu(theta, phi, t, nu)

    if jetType == jet.Spherical:
        checkCocoonArgs(**argsDict)
    else:
        checkJetArgs(**argsDict)

    # arguments are good, full steam ahead!

    if 'z' in argsDict.keys():
        z = argsDict.pop('z')
    else:
        z = 0.0

    tz = t / (1+z)
    nuz = nu * (1+z)

    # Default spreading method
    if 'spread' in argsDict:
        if argsDict['spread'] is True:
            if jetType == -2 and 'thetaCoreGlobal' in argsDict:
                argsDict['spread'] = 8
            else:
                argsDict['spread'] = 7

    Inu = np.empty(theta.shape)
    Inu.flat[:] = jet.intensity(theta.flat, phi.flat, tz.flat, nuz.flat,
                                **argsDict)

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


def checkJetArgs(**argsDict):

    for _, x in argsDict.items():
        if not np.isfinite(x):
            raise ValueError("All parameters must be finite")


    jetType = argsDict['jetType']
    specType = argsDict['specType']

    theta_obs = argsDict['thetaObs']
    E0 = argsDict['E0']
    theta_c = argsDict['thetaCore']
    n0 = argsDict['n0']
    p = argsDict['p']
    epse = argsDict['epsilon_e']
    epsB = argsDict['epsilon_B']
    xiN = argsDict['ksiN']
    dL = argsDict['dL']

    # More-or-less universal bounds
    if theta_obs < 0.0 or theta_obs > 0.5*np.pi:
        raise ValueError("theta_obs must be in [0.0, pi/2]")
    if E0 <= 0.0:
        raise ValueError("E0 must be positive")
    if theta_c <= 0.0 or theta_c > 0.5*np.pi:
        raise ValueError("theta_c must be in (0.0, pi/2]")
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

    # Bounds on optional parameters

    if jetType != jet.TopHat:
        if 'thetaWing' not in argsDict:
            raise KeyError('This jet type requires thetaWing')
        else:
            theta_w = argsDict['thetaWing']
            if (theta_w <= 0.0 or theta_w > 0.5*np.pi):
                raise ValueError("thetaWing must be in (0.0, pi/2]")

    if jetType == jet.PowerLaw or jetType == jet.PowerLawCore:
        if 'b' not in argsDict:
            raise KeyError('This jet type requires b')
        else:
            b = argsDict['b']
            if (b <= 0.0):
                raise ValueError("b must be positive")

    # Energy Injection
    if 'L0' in argsDict:
        L0 = argsDict['L0']
        if L0 < 0.0:
            raise ValueError("L0 must be non-negative")
    if 'ts' in argsDict:
        ts = argsDict['ts']
        if ts < 0.0:
            raise ValueError("ts must be non-negative")

    # Additional Luminosity
    if 'LR' in argsDict:
        LR = argsDict['LR']
        if LR < 0.0:
            raise ValueError("LR must be non-negative")
    if 'LO' in argsDict:
        LO = argsDict['LO']
        if LO < 0.0:
            raise ValueError("LO must be non-negative")
    if 'LX' in argsDict:
        LX = argsDict['LX']
        if LX < 0.0:
            raise ValueError("LX must be non-negative")

    if 'z' in argsDict:
        if argsDict['z'] < 0.0:
            raise ValueError("z must be non-negative")

    # Model Specific bounds
    if jetType == jet.Cone and argsDict['thetaCore'] > argsDict['thetaWing']:
        raise ValueError("thetaWing must be larger than thetaCore"
                         "for cone model")

    return


def checkCocoonArgs(**argsDict):

    for _, x in args.items():
        if not np.isfinite(x):
            raise ValueError("All parameters must be finite")

    u_max = argsDict['uMax']
    u_min = args['uMin']
    Er = args['Er']
    MFast = args['MFast']
    n0 = argsDict['n0']
    p = argsDict['p']
    epse = argsDict['epsilon_e']
    epsB = argsDict['epsilon_B']
    xiN = argsDict['ksiN']
    dL = argsDict['dL']

    if u_max <= 0.0:
        raise ValueError("u_max must be positive")
    if u_min <= 0.0:
        raise ValueError("u_min must be positive")
    if Er <= 0.0:
        raise ValueError("Er must be positive")
    if MFast <= 0.0:
        raise ValueError("MFast must be positive")
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

    if 'z' in argsDict:
        if argsDict['z'] < 0.0:
            raise ValueError("z must be non-negative")

    # Energy Injection
    if 'L0' in argsDict:
        L0 = argsDict['L0']
        if L0 < 0.0:
            raise ValueError("L0 must be non-negative")
    if 'ts' in argsDict:
        ts = argsDict['ts']
        if ts < 0.0:
            raise ValueError("ts must be non-negative")

    # Additional Luminosity
    if 'LR' in argsDict:
        LR = argsDict['LR']
        if LR < 0.0:
            raise ValueError("LR must be non-negative")
    if 'LO' in argsDict:
        LO = argsDict['LO']
        if LO < 0.0:
            raise ValueError("LO must be non-negative")
    if 'LX' in argsDict:
        LX = argsDict['LX']
        if LX < 0.0:
            raise ValueError("LX must be non-negative")

    return
