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
/**
 * \file
 * \brief Definition of some constants needed for the treatment of atomic and
 * molecular orbitals.
 */
#ifndef FIRE_DEAMON_STATIC_NUMBERS_H
#define FIRE_DEAMON_STATIC_NUMBERS_H

/*! the number \f$\pi\f$ (approx. 3.141592653589793) */
extern const double Pi;
/*! the number \f$ (\frac{2}{\pi})^{\frac{3}{4}} \f$ */
extern const double two_div_by_pi_to_three_fourth;
/*! the number \f$ \sqrt{2} \f$ */
extern const double sqrt2;
/*! the number \f$ \sqrt{(\frac{\pi}{2})^{\frac{3}{4}}} \f$ */
extern const double sqrt_pihalf_to_3_4;

/*! a C-type array containing the integer numbers \f$ \frac{1}{\sqrt{(2i-1)!!}} \forall
 * i \wedge i>0 \wedge i<16 \f$ and the array index is i. The name is short for "One
 * Divided By Square-root of Double Factorial of Two"
 */
extern const double odbsdfo2[];
/*! a C-type array containing the factorial of the first 11 integer numbers greater zero
 */
extern const int factorial[];

/*! a C-type array containing the integer numbers \f$ \sqrt{\frac{2i+1}{4\pi}} \forall i
 * \wedge i>0 \wedge i<16 \f$ and the array index is i
 */
extern const double sqrt_two_lplus1_div4pi[];

/*! a C-type array containing the integer numbers \f$ \frac{1}{\sqrt{i!}} \forall i
 * \wedge i>0 \wedge i<16 \f$ and the array index is i
 */
extern const double one_div_sqrt_factorial[];
/*! a C-type array containing the inverse values of \a one_div_sqrt_factorial */
extern const double sqrt_factorial[];

#endif // FIRE_DEAMON_STATIC_NUMBERS_H
