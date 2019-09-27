# Numeric GRB Afterglow models

A Python 3 module to calculate GRB afterglow light curves and spectra. Details of the methods can be found in [Ryan et al 2019](https://arxiv.org/abs/1909.11691). Builds on [van Eerten & MacFadyen 2010](https://arxiv.org/abs/1006.5125) and [van Eerten 2018](https://arxiv.org/abs/1801.01848).  This code is under active development.

Documentation available at <https://afterglowpy.readthedocs.io/>

## Attribution

If you use this code in a publication, please refer to the package by name and cite "Ryan, G., van Eerten, H., Piro, L. and Troja, E., 2019, arXiv:1910.11691" [arXiv link](https://arxiv.org/abs/1909.11691).

## Features

_afterglowpy_ computes synchrotron emission from the forward shock of a relativistic blast wave.  It includes:
- Fully trans-relativistic shock evolution through a constant density medium.
- On-the-fly integration over the equal-observer-time slices of the shock surface.
- Approximate prescription for jet spreading.
- Arbitrary viewing angles.
- Angularly structured jets, ie. E(&theta;)
- Spherical velocity-stratified outflows, ie. E(u)
- Counter-jet emission.

It has limited support (these should be considered experimental) for:
- Initial energy injection
- Inverse comption spectra
- Spatially-resolved intensity maps
- Early coasting phase

It does not include (yet):
- External wind medium, ie. n &prop; r<sup>-2</sup>
- Synchrotron self-absorbtion
- Reverse shock emission

_afterglowpy_ has been calibrated to the BoxFit code ([van Eerten, van der Horst, & Macfadyen 2011](https://arxiv.org/abs/1110.5089), available at the [Afterglow Library](https://cosmo.nyu.edu/afterglowlibrary/boxfit2011.html)) and produces similar light curves for top hat jets (within 50% when same parameters are used) both on- and off-axis.  Its jet models by default do not include an initial coasting phase, which may effect predictions for early observations.

## Installation/Building


_afterglowpy_ is available via `pip`:
```bash
$ pip install afterglowpy
```

If you are working on a local copy of this repo and would like to install from source, you can the run the following from the top level directory of the project.
```bash
$ pip install -e .
```

## Using

**This interface will be updated to be more sensible in the VERY near future**

In your python code, import the library with `import afterglowpy as grb`.  

The main function of interest is`grb.fluxDensity(t, nu, jetType, specType, *pars, **kwargs)`.  See `tests/plotLC.py` for a simple example.

`jetType` can be -1 (top hat), 0 (Gaussian), 1 (Power Law w/ core), 2 (Gaussian w/ core), 3 (Cocoon), or 4 (Smooth Power Law).  

`specType` can be 0 (global cooling time, no inverse compton) or 1 (global cooling time, inverse compton).

For jet-like afterglows (`jetTypes` -2, -1, 0, 1, 2, and 4) `pars` has 14 positional arguments:
- `0 thetaV` viewing angle in radians
- `1 E0` on-axis isotropic equivalent energy in erg
- `2 thetaC` half-width of the jet core in radians (jetType specific)
- `3 thetaW` "wing" truncation angle of the jet, in radians
- `4 b` power for power-law structure, &theta;<sup>-b</sup>
- `5 L0` Fiducial luminosity for energy injection, in erg/s, typically 0.
- `6 q` Temporal power-law index for energy injection, typically 0.
- `7 ts` Fiducial time-scale for energy injection, in seconds, typically 0.
- `8 n0` Number density of ISM, in cm<sup>-3</sup>
- `9 p` Electron distribution power-law index (p>2)
- `10 epsilon_e` Thermal energy fraction in electrons
- `11 epsilon_B` Thermal energy fraction in magnetic field
- `12 xi_N` Fraction of electrons that get accelerated
- `13 d_L` Luminosity distance in cm

For cocoon-like afterglows (`jetType` 3) `pars` has 14 positional arguments:
- `0 umax` Initial maximum outflow 4-velocity
- `1 umin` Minium outflow 4-velocity
- `2 Ei` Fiducial energy in velocity distribution, E(>u) = E<sub>i</sub>  u<sup>-k</sup>.
- `3 k` Power-law index of energy velocity distribution  
- `4 Mej` Mass of material at `umax' in solar masses
- `5 L0` Fiducial luminosity for energy injection, in erg/s, typically 0.
- `6 q` Temporal power-law index for energy injection, typically 0.
- `7 ts` Fiducial time-scale for energy injection, in seconds, typically 0.
- `8 n0` Number density of ISM, in cm<sup>-3</sup>
- `9 p` Electron distribution power-law index (p>2)
- `10 epsilon_e` Thermal energy fraction in electrons
- `11 epsilon_B` Thermal energy fraction in magnetic field
- `12 xi_N` Fraction of electrons that get accelerated
- `13 d_L` Luminosity distance in cm

Keyword arguments are:
- `z` redshift (defaults to 0)
- `tRes` time resolution of shock-evolution scheme, number of sample points per decade in time
- `latRes` latitudinal resolution for structured jets, number of shells per `thetaC`
- `rtol` target relative tolerance of flux integration
- `spread` boolean (defaults to True), whether to allow the jet to spread.



