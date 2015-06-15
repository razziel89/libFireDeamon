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
#include <vector>
#include <skin_surface_deamon.h>

int main () {
    double shrink_factor = 0.99;
    int nr_atoms = 2;
    int nr_refinements = 1;
    std::vector<double> coord_radii_vec;
    coord_radii_vec.reserve(nr_atoms*4);

    coord_radii_vec.push_back(0.0);
    coord_radii_vec.push_back(0.0);
    coord_radii_vec.push_back(0.0);
    coord_radii_vec.push_back(1.25);

    coord_radii_vec.push_back(1.0);
    coord_radii_vec.push_back(1.0);
    coord_radii_vec.push_back(1.0);
    coord_radii_vec.push_back(1.25);

    std::vector<int> ivec;
    std::vector<double> dvec;
    std::vector<int> length;

    make_skin_surface(shrink_factor, nr_atoms, coord_radii_vec, &ivec, &dvec, &length, nr_refinements);

    return 0;
}
