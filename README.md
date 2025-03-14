# Numeric GRB Afterglow models

A Python 3 module to calculate GRB afterglow light curves and spectra. Details of the methods can be found in [Ryan et al 2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...896..166R/abstract) and [Ryan et al 2024](https://ui.adsabs.harvard.edu/abs/2024ApJ...975..131R/abstract).  Builds on [van Eerten & MacFadyen 2010](https://arxiv.org/abs/1006.5125) and [van Eerten 2018](https://arxiv.org/abs/1801.01848).  This code is under active development.

Documentation available at <https://afterglowpy.readthedocs.io/>

## Attribution

If you use this code in a publication, please refer to the package by name and cite "Ryan, G., van Eerten, H., Piro, L. and Troja, E., Astrophysical Journal *896*, 166 (2020)" [ADS link](https://ui.adsabs.harvard.edu/abs/2020ApJ...896..166R/abstract).  Upgrades including centroid motion, size, and the deep Newtonian phase are presented in "Ryan, G., van Eerten, H., Troja, E., Piro, L., O'Connor, B., and Ricci, R., Astrophysical Journal *975*, 131 (2024)" [ADS link](https://ui.adsabs.harvard.edu/abs/2024ApJ...975..131R/abstract).

## Acknowledgements

This work is funded in part by the European Union’s Horizon 2020 Programme under the AHEAD2020 project (grant agreement n. 871158).

## Features

_afterglowpy_ computes synchrotron emission from the forward shock of a relativistic blast wave.  It includes:
- Fully trans-relativistic shock evolution through a constant density medium.
- On-the-fly integration over the equal-observer-time slices of the shock surface.
- Approximate prescription for jet spreading.
- Arbitrary viewing angles.
- Angularly structured jets, ie. E(&theta;)
- Spherical velocity-stratified outflows, ie. E(u)
- Counter-jet emission.
- Deep Newtonian emission.
- Image moments suitable for astrometry: centroid position and image size.

It has limited support (these should be considered experimental) for:
- Initial energy injection
- Inverse comption spectra
- Early coasting phase

It does not include (yet):
- External wind medium, ie. n &prop; r<sup>-2</sup>
- Synchrotron self-absorbtion
- Reverse shock emission

_afterglowpy_ has been calibrated to the BoxFit code ([van Eerten, van der Horst, & Macfadyen 2011](https://arxiv.org/abs/1110.5089), available at the [Afterglow Library](https://cosmo.nyu.edu/afterglowlibrary/boxfit2011.html)) and produces similar light curves for top hat jets (within 50% when same parameters are used) both on- and off-axis.  Its jet models by default do not include an initial coasting phase, which may effect predictions for early observations.

## Changelog

### New in v0.8.1
- Numpy 2.0 compatibility
- `ignoreBounds` Boolean keyword argument to ignore built-in bounds checking on parameters.

### New in v0.8.0
- Image size and position via the `moment` keyword.
- Deep Newtonian spectral evolution at late times via `specType=grb.jet.DeepNewtonian`

## Installation/Building

_afterglowpy_ is available via `pip`:
```bash
$ pip install afterglowpy
```

_afterglowpy_ is compatible with Numpy v1 and v2, Python 3.8+, and runs on MacOS, Linux, and Windows.

If you are working on a local copy of this repo and would like to install from source, you can the run the following from the top level directory of the project.
```bash
$ pip install -e .
```

## Using

In your python code, import the library with `import afterglowpy as grb`.  

The main function of interest is`grb.fluxDensity(t, nu, **kwargs)`.  See `examples/plotLightCurve.py` for a simple example.

For jet-like afterglows there are up to 13 required keyword arguments:

- `jetType` an integer code setting the jet structure. It can be `grb.jet.TopHat`, `grb.jet.Gaussian`, `grb.jet.PowerLawCore`, `grb.jet.GaussianCore`, `grb.jet.Spherical`, or `grb.jet.PowerLaw`.  
- `specType` an integer code specifying flags for the emissivity function and spectrum. Can be `grb.jet.SimpleSpec` (basic spectrum with &nu;<sub>m</sub> and &nu;<sub>c</sub>), `grb.jet.DeepNewtonian`, `grb.jet.EpsEBar` to interpret `epsilon_e` as &epsilon;&#773;<sub>e</sub> = &epsilon;<sub>e</sub>(p-2)/(p-1), `grb.jet.ICCooling` (simple inverse Compton effects on the cooling frequency, experimental). Multiple options can be combined with the `|` operator.
- `thetaObs` viewing angle in radians
- `E0` on-axis isotropic equivalent energy in erg
- `thetaCore` half-width of the jet core in radians (jetType specific)
- `thetaWing` "wing" truncation angle of the jet, in radians
- `b` power for power-law structure, &theta;<sup>-b</sup>
- `n0` Number density of ISM, in cm<sup>-3</sup>
- `p` Electron distribution power-law index (p>2)
- `epsilon_e` Thermal energy fraction in electrons
- `epsilon_B` Thermal energy fraction in magnetic field
- `xi_N` Fraction of electrons that get accelerated
- `d_L` Luminosity distance in cm

Optional keyword arguments for all models are:
- `z` redshift (defaults to 0)
- `spread` boolean (defaults to True), whether to allow the jet to spread.
- `counterjet` boolean (defaults to False), whether to include the counterjet
- `moment` array (integer dtype, same shape as t and nu) which sky moment to compute.
- `ignoreBounds` boolean (defaults to False), whether to ignore the built in paramter bounds checking.
- `L0` Fiducial luminosity for energy injection, in erg/s, default 0.0.
- `q` Temporal power-law index for energy injection, default 0.0.
- `ts` Fiducial time-scale for energy injection, in seconds, default 0.
- `tRes` time resolution of shock-evolution scheme, number of sample points per decade in time
- `latRes` latitudinal resolution for structured jets, number of shells per `thetaC`
- `rtol` target relative tolerance of flux integration



