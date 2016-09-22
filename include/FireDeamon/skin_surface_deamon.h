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
#ifndef MAKE_SKIN_SURFACE_HEADER_H
#define MAKE_SKIN_SURFACE_HEADER_H

#include <vector>

void make_skin_surface(double shrink_factor, std::vector<double> coord_radii_vec, std::vector<int> *ivec, std::vector<double> *dvec, std::vector<double> *nvec, std::vector<int> *length, int nr_refinements);

#endif //MAKE_SKIN_SURFACE_HEADER_H