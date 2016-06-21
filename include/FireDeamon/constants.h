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
#ifndef FIRE_DEAMON_STATIC_NUMBERS_H
#define FIRE_DEAMON_STATIC_NUMBERS_H

extern const double Pi;
extern const double two_div_by_pi_to_three_fourth;
extern const double sqrt2;
extern const double sqrt_pihalf_to_3_4;

//this name is a short form of: one_div_by_sqrt_double_factorial_of_2aminus1
// which is 1.0/((2*a-1)!!)
extern const double odbsdfo2[];
//the factorial of the index
extern const int factorial[];
//two times the index plus 1 divided by 4*pi
extern const double sqrt_two_lplus1_div4pi[];

extern const double one_div_sqrt_factorial[];
extern const double sqrt_factorial[] ;

#endif //FIRE_DEAMON_STATIC_NUMBERS_H
