# Semi-analytic GRB Afterglow models

A Python module to calculate GRB afterglow light curves and spectra. Makes liberal use of [van Eerten & MacFadyen 2010](https://arxiv.org/abs/1006.5125) and [van Eerten 2018](https://arxiv.org/abs/1801.01848).

## Installation/Building

To build:

```bash
$ python setup.py build
```

To install for development

```bash
$ python setup.py develop
```

To install

```bash
$ python setup.py install
```

## Using

*This interface will be updated to be more sensible in the VERY near future*

In your python code, import the library with `import grbpy as grb`.  

The main function of interest is`grb.fluxDensity(t, nu, jetType, specType, *pars, **kwargs)`

`jetType` can be -1 (top hat), 0 (Gaussian), 1 (Power Law w/ core), 2 (Gaussian w/ core), 3 (Cocoon), or 4 (Smooth Power Law).  

`specType` can be 0 (global cooling time, no inverse compton) or 1 (global cooling time, inverse compton).

For jet-like afterglows (`jetTypes` -2, -1, 0, 1, 2, and 4) `pars` has 13 positional arguments:
- `0 thetaV` viewing angle in radians
- `1 E0` on-axis isotropic equivalent energy in erg
- `2 thetaC` half-width of the jet core in radians (jetType specific)
- `3 thetaW` "wing" truncation angle of the jet, in radians
- `4 L0` Fiducial luminosity for energy injection, in erg/s, typically 0.
- `5 q` Temporal power-law index for energy injection, typically 0.
- `6 ts` Fiducial time-scale for energy injection, in seconds, typically 0.
- `7 n0` Number density of ISM, in cm^{-3}
- `8 p` Electron distribution power-law index
- `9 epsilon_e` Thermal energy fraction in electrons
- `10 epsilon_B` Thermal energy fraction in magnetic field
- `11 xi_N` Fraction of electrons that get accelerated
- `12 d_L` Luminosity distance in cm

For cocoon-like afterglows (`jetType` 3) `pars` has 14 positional arguments:
- `0 umax` Initial maximum outflow 4-velocity
- `1 umin` Minium outflow 4-velocity
- `2 Ei` Fiducial energy in velocity distribution, E(>u) = Ei * u^{-k}.
- `3 k` Power-law index of energy velocity distribution  
- `4 Mej` Mass of material at `umax' in solar masses
- `5 L0` Fiducial luminosity for energy injection, in erg/s, typically 0.
- `6 q` Temporal power-law index for energy injection, typically 0.
- `7 ts` Fiducial time-scale for energy injection, in seconds, typically 0.
- `8 n0` Number density of ISM, in cm^{-3}
- `9 p` Electron distribution power-law index
- `10 epsilon_e` Thermal energy fraction in electrons
- `11 epsilon_B` Thermal energy fraction in magnetic field
- `12 xi_N` Fraction of electrons that get accelerated
- `13 d_L` Luminosity distance in cm

Keyword arguments are:
- `z` redshift
- `tRes` time resolution of shock-evolution scheme, number of sample points per decade in time
- `latRes` latitudinal resolution for structured jets, number of shells per `thetaC`
- `rtol` target relative tolerance of flux integration



