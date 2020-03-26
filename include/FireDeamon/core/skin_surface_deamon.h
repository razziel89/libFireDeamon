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
 * \brief Create a skin surface around a set of spheres.
 */
#ifndef MAKE_SKIN_SURFACE_HEADER_H
#define MAKE_SKIN_SURFACE_HEADER_H

#include <vector>

/**
 * \brief Create a skin surface of arbitrary high discretization around a set of
 * spheres.
 *
 * A definition for skin surfaces can be found here:
 * http://doc.cgal.org/latest/Skin_surface_3/index.html You can imagine a skin surface
 * as a rubber skin contracting around a set of spheres. The degree of contraction can
 * be specified to get a sharper or smoother approximation of the spheres. First, a very
 * weakly discretized surface is generated (a sphere roughly approximated by an
 * octaeder), which can then be further refined by adding a further point in the middle
 * of every edge (for each refinement step). Increasing the number of refinement steps
 * by one more than doubles the memory requirement.
 *
 * \bug crashes if \a shrink_factor is \f$\le 0\f$ or \f$\ge 1\f$
 *
 * \bug if \a nr_refinements is large (\f$\ge 4\f$ for a system with 8GB RAM), the
 * isosurface cannot be kept in memory but no error is thrown.
 */
void make_skin_surface(
    /*! \brief double - the shrink factor that defined how "tight" the skin surface
     * shall be A value closer to 1 causes a more accurate reproduction of the union of
     * the spheres.
     */
    double shrink_factor,
    /*! \brief std::vector<double> - a flat vector containing the coordinates and radii
     * For each sphere in the set, this vector contains the three Cartesian coordinates
     * of its center followed by the radius. That means this vector has a length of 4
     * times the number of spheres in the set.
     */
    std::vector<double> coord_radii_vec,
    /*! pointer to std::vector<int> - this flat vector will be filled with triples of
     * indices that specify the facets of the skin surface
     */
    std::vector<int> *ivec,
    /*! pointer to std::vector<double> - this flat vector will be filled with triples of
     * values specifying the Cartesian coordinates of the vertices of the skin surface
     */
    std::vector<double> *dvec,
    /*! pointer to std::vector<double> - this flat vector will be filled with triples of
     * values that specify the normal vectors associated with each vertex
     */
    std::vector<double> *nvec,
    /*! pointer to std::vector<int> - this flat vector will contain the number of
     * vertices and the number of facets, in that order
     */
    std::vector<int> *length,
    /*! int - the number of refinement steps to perform */
    int nr_refinements);

#endif // MAKE_SKIN_SURFACE_HEADER_H
