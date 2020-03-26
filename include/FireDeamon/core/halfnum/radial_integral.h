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
 * \brief Contains a class that allows for computing radial integrals that appear in
 * pseudopotential integrals.
 *
 * Please see the documentation for angular_integral.h for further details about the
 * maths involved. These integrals can be written as \f$T_{N}^{\lambda}\f$.
 */
#ifndef HALFNUM_RADIAL_INTEGRAL_H
#define HALFNUM_RADIAL_INTEGRAL_H

/**
 * \brief A class that allows for computing radial integrals that appear in
 * pseudopotential integrals.
 *
 * Please see the documentation for angular_integral.h and the class AngInt for further
 * details about the maths involved. These integrals are used to compute the
 * electrostatic potential at arbitrary points in space due to molecular orbitals. The
 * integrals are computed for the products of two primitive Cartesian Gaussian
 * functions.  The integrals can be written as \f$T_{N}^{\lambda}\f$.
 *
 * The integral is computed in a coordinate system that is centered at the position at
 * which the potential shall be computed. First, the integration is initialized using
 * the exponential factor \a eta and the center of the combined Gaussian \a P and a lot
 * of helper variables are initialized that allow for fast and numerically stable
 * computation of the radial integral.
 */
class RadInt {
private:
  /*! double - the norm of the vector that is the center of the Gaussian function that
   * is the product of the two original primitive Cartesian Gaussian functions.
   */
  double P;
  /*! double - \f$\eta\f$ the sum of the exponential factors (i.e., \f$\alpha\f$ in
   * \f$\mathrm{e}^{\alpha \vec{r}}\f$) of the two original primitive Cartesian
   * Gaussian functions
   */
  double eta;
  /*! \f$\eta\cot P^2\f$ */
  double etaPP;
  /*! \f$P^2\f$ */
  double PP;
  /*! \f$\eta^2\cdot P^2\f$ */
  double etaPetaP;
  /*! \f$\sqrt{\pi}\cdot\mathrm{erf}(P\cdot\sqrt{eta})\f$ */
  double erfetaP;
  /*! \f$\mathrm{e}^{-\eta\cdot P^2}\f$ */
  double expetaPP;
  /*! \f$\frac{1}{ \eta }\f$ */
  double _eta;
  /*! \f$\frac{1}{ P }\f$ */
  double _P;
  /*! \f$\frac{1}{ \sqrt{\eta} }\f$ */
  double _sqrteta;
  /*! \f$\frac{1}{ \eta^2 }\f$ */
  double _etaeta;
  /*! \f$\frac{1}{ \eta^{\frac{3}{2}} }\f$ */
  double _eta3half;
  /*! \f$\frac{1}{ \eta\cdot P }\f$ */
  double _etaP;
  /*! \f$\frac{1}{ \eta\cdot P^2 }\f$ */
  double _etaPP;
  /*! \f$\frac{1}{ \eta^2\cdot P^2 }\f$ */
  double _etaPetaP;

public:
  /**
   * \brief Initialization function for the radial integration
   * \param eta double - the exponential factor of the combined Gaussian (i.e.,
   * sum of the original ones)
   * \param P double - norm of the vector of the center of the combined Gaussian
   * function
   */
  void Init(double eta, double P);
  /**
   * \brief Compute the radial integral
   * \param N int - parameter N of the radial integral
   * \param lambda int - parameter \f$\lambda\f$ of the radial integral
   * \return the integral value
   */
  double GetRadInt(int N, int lambda);
};
#endif // HALFNUM_RADIAL_INTEGRAL_H
