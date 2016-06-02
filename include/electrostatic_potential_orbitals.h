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
#ifndef ELECTROSTATIC_POTENTIAL_ORBITALS_H
#define ELECTROSTATIC_POTENTIAL_ORBITALS_H

#include <vector>

void electrostatic_potential_orbitals(bool progress_reports, int num_primitives, std::vector<double> prim_centers, std::vector<double> prim_exponents, std::vector<double> prim_coefficients, std::vector<int> prim_angular, std::vector<double> potential_grid, std::vector<double>  P_matrix, std::vector<int> screen, std::vector<double> *potential, double cutoff);

#endif //ELECTROSTATIC_POTENTIAL_ORBITALS_H
