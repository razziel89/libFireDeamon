.. FireDeamon documentation master file, created by
   sphinx-quickstart on Thu Mar 26 12:47:07 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FireDeamon's documentation!
=========================================

.. toctree::
    :maxdepth: 2
    :caption: Contents:

Find the C++ documentation at `CPP`_!

.. _CPP: ./cpp/index.html


Introduction
============

What you are currently viewing contains the documentation for FireDeamon, a
Python package with a C++-backend written to perform some tasks related to what
I did during my time as a PhD student that will be detailed in this
documentation. The entire project is licensed under the GNU General Public
License v3.

The package FireDeamon contains functionality that I think is useful for
people working in physical chemistry or quantum chemistry/physics, who perform
quantum chemical calculations and evaluate them afterwards (or use them in any
other way). It consists of functionality that I could not find anywhere at all
or not anywhere I could just use it (like when it's in proprietary
software).

The functionality includes:
    * a generic way to compute values defined on an arbitrary grid:
        * computed from values defined on an (not necessarily identical)
          arbitrary grid
        * realized via variadic templates
        * supports progress reports during the computation
    * finding local minima in volumetric data on arbitrary grids
    * compute the following chemical/physical quantities:
        * electron densities (from atomic basis sets)
        * electrostatic potentials from
            * clouds of point charges
            * atomic basis sets
    * compute isosurfaces through volumetric data sets:
        * arbitrarily well discretized
        * only regular grids supported
    * compute skin-surfaces around a set of spheres:
        * arbitrarily well discretized
        * arbitrary radii supported
    * interpolate quantities on arbitrary grids using:
        * nearest-neighbour interpolation
        * interpolation using inverse-distance weighting:
    * compute overlaps of atomic orbitals

The package FireDeamon has been designed to be mainly used from Python via
the provided language bindings. Many of the C++ functions are not that easy to
use (i.e., their input is not that easily prepared in the correct format) and
some sanity checks are missing. In contrast to that, the high-level Python
wrapper functions perform many sanity checks and the input is more easily
prepared properly. Thus, the Python bindings are the only supported way of using
this package.


Prerequisites
=============

You need the following dependencies to use this package:
    * a C++ compiler that supports the C++14 standard
    * GNU make
    * Python:
        * `Anaconda`_-based virtual environments highly recommended
        * including the pip package manager
        * including the packages setuptools and wheel

To perform any computation involving surfaces, you also need:
    * CGAL (Computational Geometry Algorithms Library)
    * the Boost C++ libraries


Installation
============

If you have all the aforementioned dependencies installed and you environment
activated, simply run `pip install FireDeamon`.

.. _Anaconda: https://www.anaconda.com/


The Python Package
==================

.. automodule:: FireDeamon
    :members:
    :undoc-members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
