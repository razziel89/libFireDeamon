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
#include <electrostatic_potential.h>

//#include <cstdlib>
int main () {
    int num_points = 4;
    int num_charges = 1;
    int num_threads = 3;
    
    std::vector<double> points;
    std::vector<double> charges_coordinates;
    std::vector<double> potential;
    points.reserve(3*num_points);
    potential.reserve(num_points);
    charges_coordinates.reserve(4*num_charges);

    //first 3: coordinates
    charges_coordinates.push_back( 2.0);
    charges_coordinates.push_back(-2.0);
    charges_coordinates.push_back( 3.0);
    //last: charge
    charges_coordinates.push_back(1.0);

    points.push_back(1.0);
    points.push_back(1.0);
    points.push_back(1.0);

    points.push_back(1.0);
    points.push_back(0.0);
    points.push_back(1.0);

    points.push_back(1.0);
    points.push_back(0.0);
    points.push_back(0.0);

    points.push_back(0.0);
    points.push_back(0.0);
    points.push_back(1.0);

    electrostatic_potential(num_points, num_charges, points, charges_coordinates, &potential);

    return 0;
}
