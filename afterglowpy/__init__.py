#!/usr/bin/env python
"""
===========
afterglowpy
===========

Top-level module for *afterglowpy*, computes light curves and spectra of
Gamma-ray burst (GRB) afterglows.

This module provides the primary user-facing functions provided by afterglowpy,
useful physical constants, and conversion factors. Internal computations are
performed in the sub-modules.

Functions
---------

These are the primary interface into *afterglowpy*. If you just want to compute
light curves, spectra, and intensity maps, and are unconcerned with the
underlying algorithm, these are all you need.

=============== ============================================================
*fluxDensity*   Compute the flux density F_nu of a GRB afterglow.
*intensity*     Compute the specific intensity I_nu of a GRB afterglow.
=============== ============================================================

Sub-Modules
-----------

These sub-modules perform the internal computations involved with calculating
synchrotron emission from a blast wave: calculating the evolution of a blast
wave with time, computing the synchrotron emissivity, integrating over
equal-observer-time hypersurfaces, and constructing structured jets.

=============== ============================================================
*shock*         Routines for computing evolution of a relativistic shock
*jet*           Routines for computing synchrotron emission from a jet
*cocoon*        Routines for computing synctrotron emission from a spherical
                shell.
=============== ============================================================

"""
from .version import __version__
from . import shock
from . import cocoon
from . import jet
from .flux import fluxDensity, intensity
from .cocoon import (Hz2eV, Msun, c, cgs2mJy, day2sec, eV2Hz, ee, h, hbar,
                     mJy2cgs, me, mp, parsec, sec2day, sigmaT)

__all__ = ['__version__',
           'shock', 'cocoon', 'jet', 'fluxDensity', 'intensity',
           'Hz2eV', 'Msun', 'c', 'cgs2mJy', 'day2sec', 'eV2Hz', 'ee', 'h',
           'hbar', 'mJy2cgs',
           'me', 'mp', 'parsec', 'sec2day', 'sigmaT']
