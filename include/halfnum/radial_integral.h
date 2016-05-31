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
#ifndef HALFNUM_RADIAL_INTEGRAL_H
#define HALFNUM_RADIAL_INTEGRAL_H

#define MINP 1 //do not set this to anything else than 1
#define NRSTEPS 16 //this equals 65535 quadrature points
#define MAXNP1 8 //the actual max N is this minus 1 but N==0 will never be obtained
#define RADINTTOLERANCE 1.0e-10 //an integral will only be computed if after MINNRINTS steps it is above this threshold
#define MINNRINTS 3 //the minimum number of integration steps to be performed
class RadInt{
    private:
        double *m_abscissas, *m_weights;
        double *m_r_to_N;
        size_t m_nr_elements;
        double m_T[NRSTEPS];
        int m_current_T;
        unsigned int m_p;
        unsigned int m_lambda;
        double *m_absit;
        double *m_weiit;
        double *m_radit;
        double m_P, m_eta;
        double m_epsilon;

        void FillAbscissasWeightsDifferentialsRadii();
        double* GetRStart(unsigned int N) const;
        void NewQuadrature(double eta, double P, unsigned int N, unsigned int lambda);
        bool CheckConvergence();
        void NextT();

    public:
        RadInt(double epsilon);
        ~RadInt();
        double GetInt(double eta, double P, unsigned int N, unsigned int lambda);
};
#endif //HALFNUM_RADIAL_INTEGRAL_H
