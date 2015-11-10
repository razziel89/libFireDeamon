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
#ifndef ELECTRON_DENSITY_H 
#define ELECTRON_DENSITY_H

#include <vector>

void electron_density(bool progress_reports, int num_gridpoints, std::vector<double> prim_centers, std::vector<double> prim_exponents, std::vector<double> prim_coefficients, std::vector<double> prim_angular, std::vector<double> density_grid, std::vector<double>  mo_coefficients, std::vector<double> *density);

#endif //ELECTRON_DENSITY_H
