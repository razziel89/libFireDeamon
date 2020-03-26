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
 * \brief Compute the electrostatic potential due to molecular orbitals.
 *
 * These are defined as a linear combination of primitive Cartesian Gaussian
 * functions.
 */
#ifndef ELECTROSTATIC_POTENTIAL_ORBITALS_H
#define ELECTROSTATIC_POTENTIAL_ORBITALS_H

#include <vector>

/** \brief Compute the electrostatic potential due to molecular orbitals.
 *
 * Some matrices (P and S matrices) are usually computed on the level of contracted
 * Cartesian Gaussian functions. However, this functions needs them <em>spread onto the
 * primitives</em>, which means nothing more that, if a contracted function has j
 * primitives, the value has to be duplicated j times in direct succession.
 */
void electrostatic_potential_orbitals(
    /*! whether or not to print progress reports during the computation */
    bool progress_reports,
    /*! int - the number of primitive functions making up the basis */
    int num_primitives,
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
    /*! std::vector<double> - a flat vector containing the Cartesian coordinates of the
     * points at which to compute the potential
     */
    std::vector<double> potential_grid,
    /*! std::vector<double> - a flat vector containing the first order density matrix.
     * This matrix has to be <em>spread onto the primitives</em>
     */
    std::vector<double> P_matrix,
    /*! std::vector<double> - a flat vector containing the overlap matrix of the
     * contracted Cartesian Gaussian functions. This matrix has to be <em>spread onto
     * the primitives</em>
     */
    std::vector<double> S_matrix,
    /*! pointer to std::vector<double> - this vector will hold the computed potential in
     * the same order as the coordinates were defined in \a potential_grid
     */
    std::vector<double> *potential);

#endif // ELECTROSTATIC_POTENTIAL_ORBITALS_H
