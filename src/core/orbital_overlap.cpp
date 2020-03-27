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
#include <FireDeamon/core/orbital_overlap.h>
#include <math.h>

double normalization_coefficient(double alpha, int l, int m, int n) {
  return pow(2 * alpha, 0.75) * pow(4.0 * alpha, 0.5 * (l + m + n)) * odbsdfo2[l] *
         odbsdfo2[m] * odbsdfo2[n];
}

inline double power(double base, unsigned int exp) {
  double result = 1.0;
  while (exp) {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
}

inline int df(int i) {
  int result = 1;
  while (i > 1) {
    result *= i;
    i -= 2;
  }
  return result;
}

inline int binomial(int n, int k) {
  return (factorial[n]) / (factorial[n - k] * factorial[k]);
}

double Sxyz(int a, int b, double diffA, double diffB, double gamma) {
  double result = 0.0;
  for (int i = 0; i <= a; ++i) {
    for (int j = 0; j <= b; ++j) {
      if ((i + j) % 2 == 0) {
        int factor = binomial(a, i) * binomial(b, j) * df(i + j - 1);
        double Apow = power(diffA, a - i);
        double Bpow = power(diffB, b - j);
        double gamma_pow = pow(2 * gamma, 0.5 * (i + j));
        result += factor * Apow * Bpow / gamma_pow;
      }
    }
  }
  return result / sqrt(gamma);
}
