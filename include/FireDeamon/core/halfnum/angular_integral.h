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
 * \brief Contains classes that help in computing angular integrals that appear in
 * pseudopotential integrals.
 *
 * The algorithm that performs these computations is based on the following paper:
 * Flores-Moreno, R., Alvarez-Mendez, R. J., Vela, A. and Köster, A. M. (2006),
 * Half-numerical evaluation of pseudopotential integrals. J. Comput. Chem., 27:
 * 1009–1019. doi:10.1002/jcc.20410
 */
#ifndef HALFNUM_ANGULAR_INTEGRALS_H
#define HALFNUM_ANGULAR_INTEGRALS_H

#define LMAXP1 6
/**
 * \brief Class that helps computing angular integrals that appear in pseudopotential
 * integrals.
 *
 * The efficiency from this algorithm stems from the fact that all the angular integrals
 * can be precomputed and then only have to be taken from the appropriate place. This
 * class computes the integrals upon creation and provides a function to then access the
 * data.
 */
class AngInt {
private:
  /*! double C-array - contains the pretabulated values */
  double *m_integrals;

public:
  //! \brief Constructor (angular integrals are computed here)
  AngInt();
  /*! \brief Access the pretabulated integral values
   *
   * The integrals can be written as \f$ \Omega^{ijk}_{00,\lambda\mu} \f$ when
   * using the notation of the provided paper (DOI: 10.1002/jcc.20410). They are
   * identical to the angular integrals that appear when computing the non-local
   * part.
   *
   * \param lambda unsigned int - the \f$\lambda\f$ index (sum of angular momenta
   * of the involved basis functions \f$+1\f$)
   * \param mu int - the \f$\mu\f$ index (magnetic quantum number of the combined
   * basis function, satisfies \f$-\lambda\le\mu\le\lambda\f$
   * \param i unsigned int - first index stemming from the expansion in unitary
   * sphere polynomials
   * \param j unsigned int - second index stemming from the expansion in unitary
   * sphere polynomials
   * \param k unsigned int - third index stemming from the expansion in unitary
   * sphere polynomials
   * \return the value of the pretabulated angular integral \f$
   * \Omega^{ijk}_{00,\lambda\mu} \f$
   */
  double GetInt(unsigned int lambda, int mu, unsigned int i, unsigned int j,
                unsigned int k) const;
  //! \brief Destructor (free all memory)
  ~AngInt();
};
#endif // HALFNUM_ANGULAR_INTEGRALS_H
