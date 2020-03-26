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
 * \brief Header defining functions for searching volumetric data for local minima.
 *
 * The search for local minima is a two-step procedure:
 * -# creation of a neighbour list
 * -# comparison of each value with those of it's associated neighbours
 * .
 * This means that a point is considered to be a local minimum if and only if
 * its associated value is smaller (you can define by how much) than those of
 * its neighbours. First, you should call one of the two functions
 * - make_neighbour_list_irregular and
 * - make_neighbour_list_regular
 * .
 * depending on what type of grid your data are defined on. Then, pass the vector
 * containing the neighbour list to local_minima_from_neighbour_list.
 */
#ifndef H_ARBITRARY_GRID_MINIMIZATION_DEAMON_
#define H_ARBITRARY_GRID_MINIMIZATION_DEAMON_

#include <vector>

/**
 * \brief Generate a list of all neighbours of an irregular grid within the given
 * cutoff.
 *
 * \bug segfault (at least undefined behaviour) if \a max_nr_neighbours is smaller than
 * the number of possible neighbours a point might have
 */
void make_neighbour_list_irregular(
    /*! bool - whether or not to give progress reports */
    bool progress_reports,
    /*! int - the total number of points in the grid */
    int nr_gridpoints,
    /*! int - a number larger than the maximum number of points within the cutoff any
     * single point might have
     */
    int max_nr_neighbours,
    /*! int - the desired number of neighbours per point */
    int nr_neighbours,
    /*! int - the desired type of metric to compute whether or not points are
     * neighbours, possible values are: - 1: nearest neighbours - 2: Manhattan metric
     * independent for all 3 Cartesian directions - 3: Manhattan metric
     */
    int cutoff_type,
    /*! std::vector<double> - a flat list of all the point Coordinates of the grid
     * (i.e.: [x1, y1, z1, x2, y2, z2, ..., xN, yN, zN] if N == \a nr_gridpoints
     */
    std::vector<double> points,
    /*! std::vector<double> - cutoff above which points are no longer considered to be
     * neighbours. If \a cutoff_type == 1 or 3, only the first entry in distance_cutoff
     * is used. Otherwise, the first three elements are used (cutoff for x, y and z
     * direction, respectively)
     */
    std::vector<double> distance_cutoff,
    /*! pointer to std::vector<int> - this vector will be filled with the neighbour
     * list, which is a flat list containing several entries. Each entry consists of the
     * index of a point followed by the indices of its neighbours. If an entry is -1, it
     * is to be ignored.
     */
    std::vector<int> *neighbour_list,
    /*! bool - whether or not to sort each point's neighbours by their distance from it.
     * BEWARE: when set to \a false, you might not get the nearest neighbours if \a
     * max_nr_neighbours > \a nr_neighbours
     */
    bool sort_it = true);
/**
 * \brief Generate a list of all neighbours of a regular grid within the given cutoff.
 *
 * Although the parameters are called \a nr_gridpoints_x, \a nr_gridpoints_y and \a
 * nr_gridpoints_z, grids whose axes are not perpendicular to each other can also be
 * treated (by just calling the actual axes x, y and z). All explanations here, however,
 * for the sake of simplicity, assume a cubic grid
 */
void make_neighbour_list_regular(
    /*! bool - whether or not to give progress reports */
    bool progress_reports,
    /*! bool - whether or not to allow points close to the border to be possible
     * candidates for minima
     */
    bool exclude_border,
    /*! int - how many points in the first direction the regular grid has */
    int nr_gridpoints_x,
    /*! int - how many points in the second direction the regular grid has */
    int nr_gridpoints_y,
    /*! int - how many points in the third direction the regular grid has */
    int nr_gridpoints_z,
    /*! int - let \a p be the point we look at, then find all points that lie within a
     * cube whose side length is two times \a nr_neighbour_shells the grid's lattice
     * constant (e.g., 1 means all 26 points on the first enclosing cube)
     */
    int nr_neighbour_shells,
    /*! pointer to std::vector<int> - this vector will be filled with the neighbour
     * list, which is a flat list containing several entries. Each entry consists of the
     * index of a point followed by the indices of its neighbours. If an entry is -1, it
     * is to be ignored.
     */
    std::vector<int> *neighbour_list);
/**
 * \brief Extract the indices of local minimum points using a pre-computed neighbour
 * list.
 *
 * A local minimum is defined as a point whose associated value is smaller than
 * that of all surrounding points (given the degeneration cutoff).  Setting a
 * negative degeneration cutoff means that a point has to have an associated
 * value at least the absolute value of the given degeneration cutoff smaller
 * than any sourrounding point to be considered a minimum.
 */
void local_minima_from_neighbour_list(
    /*! bool - whether or not to give progress reports */
    bool progress_reports,
    /*! int - the number of neighbours each point has (used to separate entries in \a
     * neighbour_list
     */
    int nr_neighbours,
    /*! int - \a nr_values times \a nr_neighbours must be the length of \a
     * neighbour_list
     */
    int nr_values,
    /*! std::vector<int> - what make_neighbour_list_irregular or
     * make_neighbour_list_regular fill
     */
    std::vector<int> neighbour_list,
    /*! std::vector<double> - the values associated with each point on the grid. If an
     * irregular grid was used, they have to be in the same order as the points that
     * were given to make_neighbour_list_irregular.
     */
    std::vector<double> values,
    /*! pointer to std::vector<int> - this will be filled with the indices of those
     * points that are local minima
     */
    std::vector<int> *minima,
    /*! std::vector<double> - the first value will be used a s a degeneration cutoff,
     * i.e., a point's associated value has to be this much larger than that of its
     * neighbours to be considered a local minimum (can be negative)
     */
    std::vector<double> degeneration_cutoffs,
    /*! bool - whether or not to use the value in \a upper_cutoff */
    bool use_upper_cutoff = false,
    /*! bool - whether or not to use the value in \a lower_cutoff */
    bool use_lower_cutoff = false,
    /*! double - a point whose associated value is above this number can never be a
     * minimum
     */
    double upper_cutoff = 0.0,
    /*! double - a point whose associated value is below this number can never be a
     * minimum
     */
    double lower_cutoff = 0.0,
    /*! bool - whether or not to sort the resulting minima by their depth */
    int sort_it = 0,
    /*! pointer to std::vector<double> - if not NULL, fill this vector with the depth of
     * the minima (how much "lower" their values are than that of their neighbours)
     */
    std::vector<double> *depths = NULL);

#endif // H_ARBITRARY_GRID_MINIMIZATION_DEAMON_
