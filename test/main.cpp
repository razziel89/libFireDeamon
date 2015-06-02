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
#include <iostream>

void return_polyhedron(int length_vert, int length_face, d_vec_t dvec, i_vec_t ivec) {
    //output/return variables
    std::cout << length_vert << " " << length_face << std::endl;

    for (d_vec_t_it it = dvec.begin(); it!=dvec.end();) {
//        std::cout << "(" << *it++ << ",";
//        std::cout << *it++ << ",";
//        std::cout << *it++ << ")" << std::endl;
        std::cout << *it++ << " ";
        std::cout << *it++ << " ";
        std::cout << *it++ << std::endl;
    }

    for (i_vec_t_it it = ivec.begin(); it!=ivec.end();) {
        std::cout << *it++ << " ";
        std::cout << *it++ << " ";
        std::cout << *it++ << std::endl;
    }
}

int main() {

    //these will be input variables later
    double shrink_factor = 0.5;
    int nr_atoms = 4;
    d_vec_t coord_radii_vec;
    coord_radii_vec.reserve(4*nr_atoms);
    coord_radii_vec.push_back( 1.0);
    coord_radii_vec.push_back(-1.0);
    coord_radii_vec.push_back(-1.0);
    coord_radii_vec.push_back(1.25);

    coord_radii_vec.push_back( 1.0);
    coord_radii_vec.push_back( 1.0);
    coord_radii_vec.push_back( 1.0);
    coord_radii_vec.push_back(1.25);

    coord_radii_vec.push_back(-1.0);
    coord_radii_vec.push_back( 1.0);
    coord_radii_vec.push_back(-1.0);
    coord_radii_vec.push_back(1.25);

    coord_radii_vec.push_back(-1.0);
    coord_radii_vec.push_back(-1.0);
    coord_radii_vec.push_back( 1.0);
    coord_radii_vec.push_back(1.25);

    int length_vert=0;
    int length_face=0;
    d_vec_t dvec;
    i_vec_t ivec;

    //create skin surface in a nice way
    make_skin_surface(shrink_factor, nr_atoms, coord_radii_vec, dvec, ivec, &length_vert, &length_face);
    
    //output result
    return_polyhedron(length_vert, length_face, dvec, ivec);
    
    //end the function by returning
    return 0;
}
