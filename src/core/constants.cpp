/***********
This file is part of libFireDeamon.

Copyright (C) 2016 by Torsten Sachse

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
#include <FireDeamon/core/constants.h>
#include <math.h>

const double Pi = acos(-1.0L);
const double two_div_by_pi_to_three_fourth = 0.71270547035499016035;
const double sqrt2 = sqrt(2.0L);
const double sqrt_pihalf_to_3_4 = 1.4031041455342160267;

// this name is a short form of: one_div_by_sqrt_double_factorial_of_2aminus1
// which is 1.0/sqrt((2*a-1)!!)
const double odbsdfo2[] = {1.0 / sqrt(1.0),
                           1.0 / sqrt(1.0),
                           1.0 / sqrt(3.0),
                           1.0 / sqrt(15.0),
                           1.0 / sqrt(105.0),
                           1.0 / sqrt(945.0),
                           1.0 / sqrt(10395.0),
                           1.0 / sqrt(135135.0),
                           1.0 / sqrt(2027025.0),
                           1.0 / sqrt(34459425.0),
                           1.0 / sqrt(654729075.0),
                           1.0 / sqrt(13749310575.0),
                           1.0 / sqrt(316234143225.0),
                           1.0 / sqrt(7905853580625.0),
                           1.0 / sqrt(213458046676875.0),
                           1.0 / sqrt(6190283353629375.0),
                           1.0 / sqrt(191898783962510625.0),
                           1.0 / sqrt(6332659870762850625.0)};

const int factorial[] = {
    1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600};

const double sqrt_two_lplus1_div4pi[] = {
    sqrt(1 / (4.0 * Pi)),  sqrt(3 / (4.0 * Pi)),  sqrt(5 / (4.0 * Pi)),
    sqrt(7 / (4.0 * Pi)),  sqrt(9 / (4.0 * Pi)),  sqrt(11 / (4.0 * Pi)),
    sqrt(13 / (4.0 * Pi)), sqrt(15 / (4.0 * Pi)), sqrt(17 / (4.0 * Pi)),
    sqrt(19 / (4.0 * Pi)), sqrt(21 / (4.0 * Pi)), sqrt(23 / (4.0 * Pi)),
    sqrt(25 / (4.0 * Pi)), sqrt(27 / (4.0 * Pi)), sqrt(29 / (4.0 * Pi)),
    sqrt(31 / (4.0 * Pi)), sqrt(33 / (4.0 * Pi)), sqrt(35 / (4.0 * Pi)),
    sqrt(37 / (4.0 * Pi)), sqrt(39 / (4.0 * Pi))};

const double one_div_sqrt_factorial[] = {1.0 / (sqrt(1.0)),
                                         1.0 / (sqrt(1.0)),
                                         1.0 / (sqrt(2.0)),
                                         1.0 / (sqrt(6.0)),
                                         1.0 / (sqrt(24.0)),
                                         1.0 / (sqrt(120.0)),
                                         1.0 / (sqrt(720.0)),
                                         1.0 / (sqrt(5040.0)),
                                         1.0 / (sqrt(40320.0)),
                                         1.0 / (sqrt(362880.0)),
                                         1.0 / (sqrt(3628800.0)),
                                         1.0 / (sqrt(39916800.0)),
                                         1.0 / (sqrt(479001600.0))};

const double sqrt_factorial[] = {sqrt(1.0),
                                 sqrt(1.0),
                                 sqrt(2.0),
                                 sqrt(6.0),
                                 sqrt(24.0),
                                 sqrt(120.0),
                                 sqrt(720.0),
                                 sqrt(5040.0),
                                 sqrt(40320.0),
                                 sqrt(362880.0),
                                 sqrt(3628800.0),
                                 sqrt(39916800.0),
                                 sqrt(479001600.0)};
