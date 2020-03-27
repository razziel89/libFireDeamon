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
#include <FireDeamon/core/halfnum/radial_integral.h>
#include <cmath>
#include <stdexcept>

static const double sqrtPi = sqrt(acos(-1.0L));

template <typename T> inline T power(T base, unsigned int exp) {
  T result = 1.0;
  while (exp) {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
}

void RadInt::Init(double eta_in, double P_in) {
  eta = eta_in;
  P = P_in;
  etaPP = eta * P * P;
  PP = P * P;
  double sqrteta = sqrt(eta);
  erfetaP = sqrtPi * erf(sqrteta * P);
  expetaPP = exp(-etaPP);
  _eta = 1.0 / eta;
  _sqrteta = 1.0 / sqrteta;
  _P = 1.0 / P;
  _etaeta = _eta * _eta;
  _eta3half = _eta * _sqrteta;
  _etaP = _eta * _P;
  _etaPP = _etaP * _P;
  _etaPetaP = _etaP * _etaP;
}

double RadInt::GetRadInt(int N, int lambda) {
  double result = 0.0;
  switch (N) {
  case 1:
    switch (lambda) {
    case 0:
      result = 0.25 * erfetaP * _P * _eta3half;
      break;
    default:
      throw std::runtime_error("Wrong input value for lambda.");
      break;
    }
    break;
  case 2:
    switch (lambda) {
    case 0:
      result = 0.25 * sqrtPi * _eta3half;
      break;
    case 1:
      result = 0.25 * expetaPP * _etaP * _eta - 0.125 * erfetaP * _etaPetaP * _sqrteta +
               0.25 * erfetaP * _eta3half;
      break;
    default:
      throw std::runtime_error("Wrong input value for lambda.");
      break;
    }
    break;
  case 3:
    switch (lambda) {
    case 0:
      result = 0.25 * expetaPP * _etaeta + 0.125 * erfetaP * _etaP * _eta3half +
               0.25 * P * erfetaP * _eta3half;
      break;
    case 1:
      result = 0.25 * sqrtPi * P * _eta3half;
      break;
    case 2:
      result = -0.375 * expetaPP * _etaPetaP * _eta + 0.25 * expetaPP * _etaeta +
               0.1875 * erfetaP * _etaPetaP * _etaP * _sqrteta -
               0.25 * erfetaP * _etaP * _eta3half + 0.25 * P * erfetaP * _eta3half;
      break;
    default:
      throw std::runtime_error("Wrong input value for lambda.");
      break;
    }
    break;
  case 4:
    switch (lambda) {
    case 0:
      result = 0.375 * sqrtPi * _eta * _eta3half + 0.25 * PP * sqrtPi * _eta3half;
      break;
    case 1:
      result = 0.125 * expetaPP * _etaP * _etaeta + 0.25 * expetaPP * P * _etaeta -
               0.0625 * erfetaP * _etaPetaP * _eta3half +
               0.25 * erfetaP * _eta * _eta3half + 0.25 * PP * erfetaP * _eta3half;
      break;
    case 2:
      result = 0.25 * PP * sqrtPi * _eta3half;
      break;
    case 3:
      result = 0.9375 * expetaPP * _etaPetaP * _etaP * _eta -
               0.5 * expetaPP * _etaP * _etaeta + 0.25 * expetaPP * P * _etaeta -
               0.46875 * erfetaP * _etaPetaP * _etaPetaP * _sqrteta +
               0.5625 * erfetaP * _etaPetaP * _eta3half -
               0.375 * erfetaP * _eta * _eta3half + 0.25 * PP * erfetaP * _eta3half;
      break;
    default:
      throw std::runtime_error("Wrong input value for lambda.");
      break;
    }
    break;
  case 5:
    switch (lambda) {
    case 0:
      result = 0.625 * expetaPP * _etaeta * _eta + 0.25 * expetaPP * PP * _etaeta +
               0.1875 * erfetaP * _etaP * _etaeta * _sqrteta +
               0.75 * P * erfetaP * _eta * _eta3half +
               0.25 * PP * P * erfetaP * _eta3half;
      break;
    case 1:
      result =
          0.625 * P * sqrtPi * _eta * _eta3half + 0.25 * PP * P * sqrtPi * _eta3half;
      break;
    case 2:
      result = -0.1875 * expetaPP * _etaPetaP * _etaeta +
               0.25 * expetaPP * _etaeta * _eta + 0.25 * expetaPP * PP * _etaeta +
               0.09375 * erfetaP * _etaPetaP * _etaP * _eta3half -
               0.1875 * erfetaP * _etaP * _eta3half * _eta +
               0.375 * P * erfetaP * _eta * _eta3half +
               0.25 * PP * P * erfetaP * _eta3half;
      break;
    case 3:
      result = 0.25 * PP * P * sqrtPi * _eta3half;
      break;
    case 4:
      result = -3.28125 * expetaPP * _etaPetaP * _etaPetaP * _eta +
               1.5625 * expetaPP * _etaPetaP * _etaeta -
               0.625 * expetaPP * _etaeta * _eta + 0.25 * expetaPP * PP * _etaeta +
               1.640625 * erfetaP * _etaPetaP * _etaPetaP * _etaP * _sqrteta -
               1.875 * erfetaP * _etaPetaP * _etaP * _eta3half +
               1.125 * erfetaP * _etaP * _eta3half * _eta -
               0.5 * P * erfetaP * _eta * _eta3half +
               0.25 * PP * P * erfetaP * _eta3half;
      break;
    default:
      throw std::runtime_error("Wrong input value for lambda.");
      break;
    }
    break;
  default:
    throw std::runtime_error("Wrong input value for N.");
  }

  return result;
}
