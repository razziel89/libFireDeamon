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
 * \brief Functions to quickly compute normalization coefficients and overlaps of
 * Cartesian Gaussian orbitals.
 */
#ifndef ORBITAL_OVERLAP_H
#define ORBITAL_OVERLAP_H

/**
 * \brief Compute the normalization factor for a primitive Cartesian Gaussian orbital
 *
 * The orbital is of the form:
 *
 * \f$ G(\vec{r}) = (x-X_0)^l (y-Y_0)^m (z-Z_0)^n \cdot \mathrm{e}^{-\alpha
 * (\vec{r}-\vec{R}_0)^2} \f$
 *
 * with
 *
 * \f$ \vec{R}_0 = \left(\begin{array}{c}X_0 \\ Y_0 \\ Z_0 \end{array}\right)\f$
 *
 * and
 *
 * \f$ \vec{r} = \left(\begin{array}{c}x \\ y \\ z \end{array}\right)\f$
 *
 * \param alpha double - the exponential factor \f$\alpha\f$
 * \param l     int    - first angular momentum factor l
 * \param m     int    - second angular momentum factor m
 * \param n     int    - third angular momentum factor n
 * \return normalization coefficient
 */
double normalization_coefficient(double alpha, int l, int m, int n);

/**
 * \brief compute the overlap between two one-dimensional Cartesian Gaussian functions
 *
 * Such functions have the form \f$ G(x) = (x-X_0)^{a/b} \cdot \mathrm{e}^{-\alpha/\beta
 * (x-X_0)^2} \f$. Such a computation can be simplified if both Gaussians are combined
 * to one Gaussian and are regarded in a coordinate system whose origin is at the center
 * of the combined Gaussian (i.e., the product of the two original ones).
 *
 * \param a int - the Cartesian factor in front of the first Cartesian Gaussian function
 * \param a int - the Cartesian factor in front of the second Cartesian Gaussian
 * function
 * \param diffA double - the difference between the center of the combined Gaussian and
 * the center of the first Gaussian
 * \param diffB double - the difference between the center of the combined Gaussian and
 * the center of the second Gaussian
 * \param gamma double - the exponent of the combined Gaussian computes as \f$ \alpha +
 * \beta \f$
 */
double Sxyz(int a, int b, double diffA, double diffB, double gamma);

#endif // ORBITAL_OVERLAP_H
