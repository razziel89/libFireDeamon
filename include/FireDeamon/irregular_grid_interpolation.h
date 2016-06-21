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
#ifndef H_IRREGULAR_GRID_INTERPOLATION_DEAMON_
#define H_IRREGULAR_GRID_INTERPOLATION_DEAMON_

#include <vector>

void generic_interpolation(bool progress_reports, int num_interpolation_points, std::vector<double> points, std::vector<double> values, std::vector<double> interpolation_points, std::vector<double> *interpolation, int interpolation_type, int distance_exponent, int distance_function, double cutoff=-1.0);

#endif //H_IRREGULAR_GRID_INTERPOLATION_DEAMON_
