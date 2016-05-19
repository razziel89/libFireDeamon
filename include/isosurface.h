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
#ifndef MAKE_ISOSURFACE_HEADER_H
#define MAKE_ISOSURFACE_HEADER_H

#include <vector>

void make_isosurface(std::vector<double> data, std::vector<double> origin,
        std::vector<double> voxel, std::vector<int> extent, std::vector<double>
        points_inside, std::vector<double> mesh_criteria, std::vector<double>
        radii, double relative_precision, double isovalue,std::vector<int> *ivec,
        std::vector<double> *dvec, std::vector<double> *nvec, std::vector<int> *length);

#endif //MAKE_ISOSURFACE_HEADER_H
