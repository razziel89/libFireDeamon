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
 * \brief Compute the electrostatic potential due to a point cloud of charges.
 */
#ifndef ELECTROSTATIC_POTENTIAL_CHARGES_H
#define ELECTROSTATIC_POTENTIAL_CHARGES_H

#include <vector>

//! \brief Compute the electrostatic potential due to a point cloud of charges.
void electrostatic_potential(
    /*! bool - whether or not to print progress reports during the computation */
    bool progress_reports,
    /*! int - at how many points shall the potential be computed */
    int num_points,
    /*! std::vector<double> - a flat list of the Cartesian coordinates of the points at
     * which to compute the potential
     */
    std::vector<double> points,
    /*! std::vector<double> - a flat list containing the information about the point
     * cloud. Each charge in the cloud is described by four values: -# its charge -# its
     * x-coordinate -# its y-coordinate -# its z-coordinate
     */
    std::vector<double> charges_coordinates,
    /*! pointer to std::vector<double> - this vector will hold the computed potential in
     * the same order as the points were specified in \a points
     */
    std::vector<double> *potential,
    /*! double - if a charge is farther away than this from a point at which the the
     * potential is to be computed, do not consider this charge. A negative value
     * switches off this behaviour.
     */
    double cutoff);

#endif // ELECTROSTATIC_POTENTIAL_CHARGES_H
