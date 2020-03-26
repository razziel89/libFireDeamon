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
 * \brief Function to create an isosurface of arbitrary high quality through volumetric
 * data.
 *
 * The function \a make_isosurface has mainly been designed to create isosurfaces arount
 * molecules. It is fast for single molecules but might take longer for multiple
 * molecules (i.e., in the case of non-overlapping isosurfaces) and might not finish if
 * certain conditions are not met.  See bugs.
 *
 * HINT: one isosurface computation is performed for each point specified in \a
 * points_inside. So declaring only the required minimum (1 in the case of a single
 * molecule) greatly speeds up the computation.
 *
 * \bug The algorithm does not yield the correct iso surface if the points declared in
 * \a points_inside are not actually located near the isosurface (they don't have to be
 * inside, but they need to be close). This bug is no problem for molecules since its
 * atoms should lie inside the isosurface.
 *
 * \bug The algorithm does not finish if the angular bound mesh criterion (first entry
 * in \a mesh_criteria) smaller than 30.0 degrees.
 *
 * \bug The algorithm does not finish if the radii given in \a radii do not define
 * spheres that completely enclose the to-be-generated isosurfaces.
 */
#ifndef MAKE_ISOSURFACE_HEADER_H
#define MAKE_ISOSURFACE_HEADER_H

#include <vector>

void make_isosurface(
    /*! std::vector<double> - a flat list containing the volumetric data. The order for
     * the indices of the data is: z - fast, y - medium, x - slow
     */
    std::vector<double> data,
    /*! std::vector<double> - a flat list containing the origin point of the data (3
     * values)
     */
    std::vector<double> origin,
    /*! std::vector<double> - a flat list containing the lengths of the voxel sides.
     * This must contain 3 values for x, y and z directions. This means that the voxel
     * vectors need to be parallel to the 3 Cartesian axes. Of course, also non-cuboid
     * voxels can be treated after mapping them to rectangular voxels.
     */
    std::vector<double> voxel,
    /*! std::vector<int> - a flat list containing the number of points in x, y and z
     * directions.
     */
    std::vector<int> extent,
    /*! std::vector<double> - a flat list containing the Cartesian coordinates for the
     * points that lie within the isosurfaces. The length has to be divisible by 3.
     */
    std::vector<double> points_inside,
    /*! std::vector<double> - a flat list containing the three meshing criteria: -#
     * Angular bound for surface mesh generation. If <30, the algorithm is not
     * guaranteed to finish. This is the lower bound in degrees for the angles during
     * mesh generation. -# Radius bound used during mesh generation. It is an upper
     * bound on the radii of surface Delaunay balls. A surface Delaunay ball is a ball
     * circumscribing a mesh facet and centered on the surface. -# Distance bound used
     * during surface mesh generation. It is an upper bound for the distance between the
     * circumcenter of a mesh facet and the center of a surface Delaunay ball of this
     * facet.
     */
    std::vector<double> mesh_criteria,
    /*! std::vector<double> - a flat list containing radii that, together with the
     * points given in \a points_inside, define spheres that MUST completely enclose the
     * isosurface that will be generated. I recommend choosing values large enough so
     * that the entire volumetric data set is enclosed.
     */
    std::vector<double> radii,
    /*! double - precision value used to compute the isosurface (given relative to the
     * radii). A lower value results in more highly discretized isosurfaces.
     */
    double relative_precision,
    /*! double - the isovalue at which to compute the isosurface */
    double isovalue,
    /*! pointer to std::vector<int> - this flat vector will be filled with triples of
     * indices that specify the facets of the isosurface
     */
    std::vector<int> *ivec,
    /*! pointer to std::vector<double> - this flat vector will be filled with triples of
     * values specifying the Cartesian coordinates of the vertices of the isosurface
     */
    std::vector<double> *dvec,
    /*! pointer to std::vector<double> - this flat vector will be filled with triples of
     * values that specify the normal vectors associated with each vertex
     */
    std::vector<double> *nvec,
    /*! pointer to std::vector<int> - this flat vector will contain the number of
     * vertices and the number of facets, in that order
     */
    std::vector<int> *length);

#endif // MAKE_ISOSURFACE_HEADER_H
