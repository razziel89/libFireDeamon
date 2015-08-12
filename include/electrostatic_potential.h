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
#ifndef ELECTROSTATIC_POTENTIAL_H
#define ELECTROSTATIC_POTENTIAL_H

#include <vector>

void electrostatic_potential (bool progress_reports, int num_points, int num_charges, std::vector<double> points, std::vector<double> charges_coordinates, std::vector<double> *potential);

#endif //ELECTROSTATIC_POTENTIAL_H
