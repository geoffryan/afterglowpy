from . import cocoon
from . import jet

def fluxDensity(t, nu, jetType, specType, *args, **kwargs):

    if jetType == 3:
        return cocoon.fluxDensity(t, nu, jetType, specType, *args, **kwargs)
    else:
        return jet.fluxDensity(t, nu, jetType, specType, *args, **kwargs)
