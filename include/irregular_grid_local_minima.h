/***********
This file is part of libFireDeamon.

Copyright (C) 2015,2016 by Torsten Sachse

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
#ifndef H_IRREGULAR_GRID_MINIMIZATION_DEAMON_
#define H_IRREGULAR_GRID_MINIMIZATION_DEAMON_

#include <vector>

//Generate a list of all neighbours within the given cutoff.
//Make sure max_nr_neighbours is greater than the number of possible neighbours a point might have.
void make_neighbour_list(bool progress_reports, int nr_gridpoints, int max_nr_neighbours, int nr_neighbours, int cutoff_type, std::vector<double> points, std::vector<double> distance_cutoff, std::vector<int>* neighbour_list, bool sort_it=true);
//A local minimum is defined as a point whose associated value is smaller than that of
//all surrounding points (given the degeneration cutoff).
//Setting a negative degeneration cutoff means that a point has to have an associated value at least
//the absolute value of the given degeneration cutoff smaller than any sourrounding point
//to be considered a minimum
void local_minima_from_neighbour_list(bool progress_reports, int nr_neighbours, int nr_values, std::vector<int> neighbour_list, std::vector<double> values, std::vector<int>* minima, std::vector<double> degeneration_cutoffs, bool use_upper_cutoff=false, bool use_lower_cutoff=false, double upper_cutoff=0.0, double lower_cutoff=0.0, int sort_it=0, std::vector<double>* depths=NULL);

#endif //H_IRREGULAR_GRID_MINIMIZATION_DEAMON_
