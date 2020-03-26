/***********
This file is part of libFireDeamon.

Copyright (C) 2016 by Torsten Sachse

libFireDeamon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libFireDeamon is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with libFireDeamon.  If not, see <http://www.gnu.org/licenses/>.
***********/
/**
 * \file
 * \brief Routines to compute the electron density as well as the overlap between
 * Gaussian-type atomic orbitals.
 */
#ifndef ELECTRON_DENSITY_H
#define ELECTRON_DENSITY_H

#include <vector>

//! \brief Compute the electron density on an arbitrary grid caused by molecular
//! orbitals.
//!
//! Molecular orbitals are given as a linear combination of atomic orbitals and
//! occupation numbers. The basis has to be specified in terms of normalized, primitive
//! Cartesian Gaussian orbitals, which means that \a prim_centers, \a prim_exponents, \a
//! prim_coefficients and \a prim_angular have to have the exact same length
//! (considering  that each primitive has one center, exponent and coefficient, but its
//! angular momentum and center in space are each described by three values).
void electron_density(
    /*! bool - whether or not to output progress reports during the computation */
    bool progress_reports,
    /*! int - the number of points at which to compute the density */
    int num_gridpoints,
    /*! std::vector<double> - a flat list of the Cartesian coordinates of the
     * primitives' center (length==3N with N==no. of primitives)
     */
    std::vector<double> prim_centers,
    /*! std::vector<double> - a flat list of the exponential factors of the primitives
     */
    std::vector<double> prim_exponents,
    /*! std::vector<double> - a flat list of the preexponential factors of the
     * primitives
     */
    std::vector<double> prim_coefficients,
    /*! std::vector<int> - a flat list of the angular factors of the Cartesian
     * primitives (length==3N with N==no. of primitives)
     */
    std::vector<int> prim_angular,
    /*! std::vector<double> - a flat list of the Cartesian coordinates at which to
     * compute the density
     */
    std::vector<double> density_grid,
    /*! std::vector<double> - a flat list of coefficients specifying how the atomic
     * basis described with the above parameters consitutes a molecular orbital
     */
    std::vector<double> mo_coefficients,
    /*! pointer to std::vector<double> - this vector will hold the resulting density
     * values
     */
    std::vector<double> *density,
    /*! double - if the center of two primitives are farther away from each other than
     * this value, do not compute the density due to the overlap of these orbitals
     */
    double cutoff = -1.0);

//! \brief Compute the normalization coefficients for a set of primitive Cartesian
//! Gaussian functions
void normalize_gaussians(
    /*! pointer to std::vector<double> - this vector will hold the computed
     * normalization coefficients in the same order used for \a exponent and \a angular
     */
    std::vector<double> *prefactor,
    /*! std::vector<double> - a flat list of the exponential factors of the primitive
     * Cartesian Gaussian functions
     */
    std::vector<double> exponent,
    /*! std::vector<int> - a flat list of the angular factors of the Cartesian
     * primitives (length==3N with N==no. of Cartesian Gaussian functions)
     */
    std::vector<int> angular);

#endif // ELECTRON_DENSITY_H
