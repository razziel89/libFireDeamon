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
#include <algorithm>
#include <cstring>
#include <cmath>

#include <halfnum/kfunction.h>

static const double Pi     = acos(-1.0L);
static const double logof2 = log(2.0);

double power(double base, unsigned int exp){
    double result = 1.0;
    while (exp){
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}
double power(int base, unsigned int exp){
    double result = 1.0;
    while (exp){
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

#define MINP 1 //do not set this to anything else than 1
#define NRSTEPS 10 //this equals 1023 quadrature points
#define MAXNP1 8 //the actual max N is this minus 1 but N==0 will never be obtained
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

        void FillAbscissasWeightsDifferentialsRadii(){
            //BEWARE: this computes the abscissas for the successive quadratures in one array!!!
            //Therefore, consider that the weights for the first quadratures need to be adjusted to be used with
            //the last one!!!
            //BEWARE: since the original quadrature integrates over x in [-1,1] but this one integrates
            //over r in 0,inf the differential is different (dr/dx = pow(2,r-1)/Log[2]). This differential is
            //computed here as well.
            double xk, wk, arg;
            int p = MINP;
            int count = 0;
            for (int k=1; k<=p; ++k, ++count){
                arg = k*Pi/(1.0*p+1.0);
                double sinarg = sin(arg);
                xk = (p+1-2*k)/(1.0*p+1.0) +
                    (2.0/Pi)*(1.0+2.0*sinarg*sinarg/3.0)*cos(arg)*sinarg;
                m_abscissas[count] = 1.0-log(1.0-xk)/logof2;
                wk = 16.0/(3.0*(p+1))*sinarg*sinarg*sinarg*sinarg;
                m_weights[count] = wk / ((1.0-xk)*logof2);
            }
            for (int step=1; step<NRSTEPS; ++step){
                p = 2*p+1;
                for (int k=1; k<=p; k+=2, ++count){
                    arg = k*Pi/(1.0*p+1.0);
                    double sinarg = sin(arg);
                    xk = (p+1-2*k)/(1.0*p+1.0) +
                        (2.0/Pi)*(1.0+2.0*sinarg*sinarg/3.0)*cos(arg)*sinarg;
                    m_abscissas[count] = log(2.0/(1.0-xk))/logof2;
                    wk = 16.0/(3.0*(p+1))*sinarg*sinarg*sinarg*sinarg;
                    m_weights[count] = wk / ((1.0-xk)*logof2);
                }
            }
            double* r_to_N_it = m_r_to_N;
            for (unsigned int N=1; N<MAXNP1; ++N){
                for (unsigned long int rit=0; rit<m_nr_elements; ++rit, ++r_to_N_it){
                    *r_to_N_it = power(m_abscissas[rit],N);
                }
            }
        }
        double* GetRStart(unsigned int N){
            if (N==0 || N>=MAXNP1){
                throw;
            }
            return m_r_to_N+((N-1)*m_nr_elements);
        }
        void NewQuadrature(double eta, double P, unsigned int N, unsigned int lambda){
            m_absit = m_abscissas;
            m_weiit = m_weights;
            m_radit = GetRStart(N);
            m_current_T = -1;
            m_P   = P;
            m_eta = eta;
            m_p   = MINP;
            m_lambda = lambda;
            for (int i=0; i<NRSTEPS; ++i){
                m_T[i] = 0.0;
            }
        }
        bool CheckConvergence(){
            return (
                    power(m_T[m_current_T]-m_T[m_current_T-1],2)
                                    <=
                    fabs(m_T[m_current_T]-m_T[m_current_T-2])*m_epsilon
                    );
        }
        void NextT(){
            double T = 0.0;
            if (m_current_T>=0){
                T = 0.5*m_T[m_current_T];
            }
            for (unsigned int k=1; k<=m_p; k+=2, ++m_absit, ++m_weiit, ++m_radit){
                double r = *m_absit;
                T += *m_weiit * *m_radit * Kfunction(2.0*m_eta*m_P*r,m_lambda) * exp(-m_eta*(r - m_P)*(r - m_P));
            }
            m_p = 2*m_p+1;
            ++m_current_T;
            m_T[m_current_T] = T;
        }

    public:
        RadInt(double epsilon){
            m_epsilon       = epsilon;
            m_nr_elements   = power(2,NRSTEPS)-1;
            m_abscissas     = (double*) malloc(m_nr_elements*sizeof(double));
            m_weights       = (double*) malloc(m_nr_elements*sizeof(double));
            memset(m_abscissas    ,0.0,m_nr_elements*sizeof(double));
            memset(m_weights      ,0.0,m_nr_elements*sizeof(double));
            m_r_to_N = (double*) malloc((MAXNP1-1)*m_nr_elements*sizeof(double));
            memset(m_r_to_N       ,0.0,(MAXNP1-1)*m_nr_elements*sizeof(double));
            FillAbscissasWeightsDifferentialsRadii();
        }
        ~RadInt(){
            free(m_abscissas);
            free(m_weights);
            free(m_r_to_N); }
        double GetInt(double eta, double P, unsigned int N, unsigned int lambda){
            NewQuadrature(eta, P, N, lambda);
            NextT();
            NextT();
            NextT();
            while (not(CheckConvergence()) && m_current_T<NRSTEPS-1){
                NextT();
            }
            if (not(CheckConvergence())){
                throw;
            }
            return m_T[m_current_T];
        }
};
#undef MINP
#undef NRSTEPS
