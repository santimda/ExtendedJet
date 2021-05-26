#!/usr/bin/env python3
"""
===========
afterglowpy
===========

Top-level module for *afterglowpy*, computes light curves and spectra of
Gamma-ray burst (GRB) afterglows.

This module provides the primary user-facing functions provided by afterglowpy,
useful physical constants, and conversion factors. Internal computations are
performed in the submodules.

Functions
---------

These are the primary interface into *afterglowpy*. If you just want to compute
light curves, spectra, and intensity maps, and are unconcerned with the
underlying algorithm, these are all you need.

=================== =========================================================
:func:`fluxDensity` Compute the flux density F_nu of a GRB afterglow.
:func:`intensity`   Compute the specific intensity I_nu of a GRB afterglow.
=================== =========================================================

Submodules
-----------

These submodules perform the internal computations involved with calculating
synchrotron emission from a blast wave: calculating the evolution of a blast
wave with time, computing the synchrotron emissivity, integrating over
equal-observer-time hypersurfaces, and constructing structured jets.

=========================== ===================================================
:mod:`afterglowpy.shock`    Routines for computing evolution of a relativistic
                            shock
:mod:`afterglowpy.jet`      Routines for computing synchrotron emission from a
                            jet
:mod:`afterglowpy.cocoon`   Routines for computing synctrotron emission from a
                            spherical shell.
=========================== ===================================================

"""
#from .version import __version__
from . import main
radsyn_data_run = main.radsyn.radsyn_data_run

__all__ = ['radsyn_data_run']
