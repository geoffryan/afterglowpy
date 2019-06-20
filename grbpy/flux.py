from . import cocoon
from . import jet
# import time
import numpy as np


class fluxCalc(object):
    """A generic afterglow flux calculator."""

    nu_stencil = 5
    t_stencil = 5

    def __init__(self):
        return

    def Fnu(self, t, nu, params):
        return np.zeros((t+nu).shape)

    def F(self, t, nu0, nu1, params):

        shape = ((t+nu0+nu1).shape[0], self.nu_stencil)

        Fnu = np.empty(shape)
        x = np.linspace(0, 1, self.nu_stencil)
        nu = (1-x)[None, :]*nu0[:, None] + x[None, :]*nu1[:, None]
        tt = t  # should be shape of nu

        Fnu.flat[:] = self.Fnu(tt, nu, params)

        return Fnu.sum(axis=1)


def fluxDensity(t, nu, jetType, specType, *args, **kwargs):

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

    if 'z' in kwargs.keys():
        z = kwargs['z']
    else:
        z = 0.0

    tz = t / (1+z)
    nuz = nu * (1+z)

    Inu = np.empty(theta.shape)
    Inu.flat[:] = jet.intensity(theta.flat, phi.flat, tz.flat, nuz.flat,
                                jetType, specType, *args, **kwargs)

    # K-correct the intensity
    # I'm only using the flux correction here, which leaves the angular
    # part of the intensity uncorrected. This function is approximate
    # anyways, so oh well!

    Inu *= 1+z

    return Inu
