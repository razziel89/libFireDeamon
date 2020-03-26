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
 * \brief Interpolate data defined on an arbitrary grid onto another arbitrary grid.
 */
#ifndef H_IRREGULAR_GRID_INTERPOLATION_DEAMON_
#define H_IRREGULAR_GRID_INTERPOLATION_DEAMON_

#include <vector>

//! \brief Interpolate data defined on an arbitrary grid A onto another arbitrary grid
//! B.
void generic_interpolation(
    /*! bool - whether or not to print progress reports during the computation */
    bool progress_reports,
    /*! int - the number of points of grid B */
    int num_interpolation_points,
    /*! std::vector<double> - a flat list containing the Cartesian coordinats of the
     * points on grid A
     */
    std::vector<double> points,
    /*! std::vector<double> - a list containung the values associated with the points
     * whose coordinats are in \a points (i.e., those of grid A)
     */
    std::vector<double> values,
    /*! std::vector<double> - a flat list containing the Cartesian coordinats of the
     * points on grid B
     */
    std::vector<double> interpolation_points,
    /*! pointer to std::vector<double> - a list that will contain the values associated
     * with the points on grid B (i.e., the interpolation result)
     */
    std::vector<double> *interpolation,
    /*! int - specify the type of interpolation to use. 1: nearest neighbour, 2: inverse
     * distance
     */
    int interpolation_type,
    /*! int - if using inverse-distance scaling, this is the exponent of the norm */
    int distance_exponent,
    /*! int - if using inverse-distance scaling, declare the norm to use. The number 2
     * means the Eukledian norm, 3 the 3-norm, etc.
     */
    int distance_function,
    /*! double - if a point in grid A is farther away from a point in grid B than this
     * value, do not consider the value at that A-point to get the value at the B-point
     */
    double cutoff = -1.0);

#endif // H_IRREGULAR_GRID_INTERPOLATION_DEAMON_
