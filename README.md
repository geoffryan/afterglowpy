# Semi-analytic GRB Afterglow models

A Python module to interface with Hendrik van Eerten's semi-analytic afterglow model.

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

In your python code, import the library with `import grbpy as grb`.  This main functions are `grb.fluxDensity(t, nu, jetType, specType, *pars)` for emission from a jet and `grb.fluxDensityCocoon(t, nu, jetType, specType, *pars)` for quasi-spherical emission.  For an example, look at `tests/plotLC.py`.


`jetType` can be -1 (top hat), 0 (Gaussian), 1 (Power Law w/ core), 2 (Gaussian w/ core), 3 (Cocoon), or 4 (Smooth Power Law).  

`specType` can be 0 (global cooling time, no inverse compton) or 1 (global cooling time, inverse compton).

