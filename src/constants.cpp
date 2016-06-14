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
#include <constants.h>
#include <math.h>

const double Pi = acos(-1.0L);
const double two_div_by_pi_to_three_fourth = 0.71270547035499016035;
const double sqrt2 = sqrt(2.0L);
const double sqrt_pihalf_to_3_4 = 1.4031041455342160267;

//this name is a short form of: one_div_by_sqrt_double_factorial_of_2aminus1
// which is 1.0/((2*a-1)!!)
const double odbsdfo2[] = {
    1.0/sqrt(1.0),
    1.0/sqrt(1.0),
    1.0/sqrt(3.0),
    1.0/sqrt(15.0),
    1.0/sqrt(105.0),
    1.0/sqrt(945.0),
    1.0/sqrt(10395.0),
    1.0/sqrt(135135.0),
    1.0/sqrt(2027025.0),
    1.0/sqrt(34459425.0),
    1.0/sqrt(654729075.0),
    1.0/sqrt(13749310575.0),
    1.0/sqrt(316234143225.0),
    1.0/sqrt(7905853580625.0),
    1.0/sqrt(213458046676875.0),
    1.0/sqrt(6190283353629375.0),
    1.0/sqrt(191898783962510625.0),
    1.0/sqrt(6332659870762850625.0)
};

const int factorial[] = {
    1                  ,
    1                  ,
    2                  ,
    6                  ,
    24                 ,
    120                ,
    720                ,
    5040               ,
    40320              ,
    362880             ,
    3628800            ,
    39916800           ,
    479001600
};
