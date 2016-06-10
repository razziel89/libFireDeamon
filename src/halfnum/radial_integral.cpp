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
#include <stdexcept>

#include <iostream>

#include <halfnum/kfunction.h>
#include <halfnum/radial_integral.h>

#if (MINNRINTS<3)
#error "MINNRINTS must be at least 3."
#endif

#include <gsl/gsl_sf_hyperg.h>
//double gsl_sf_hyperg_1F1(double a, double b, double x);
#include <gsl/gsl_sf_gamma.h>
//double gsl_sf_gammainv(const double x);
//double gsl_sf_gamma(const double x);

#include <iostream>

static const double Pi     = acos(-1.0L);
static const double sqrtPi = sqrt(acos(-1.0L));
static const double logof2 = log(2.0);

inline double power(double base, unsigned int exp){
    double result = 1.0;
    while (exp){
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}
inline long double power(long double base, unsigned int exp){
    double result = 1.0;
    while (exp){
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}
inline double power(int base, unsigned int exp){
    double result = 1.0;
    while (exp){
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

void RadInt::FillAbscissasWeightsDifferentialsRadii(){
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
double* RadInt::GetRStart(unsigned int N) const {
    if (N==0 || N>=MAXNP1){
        if (N==0){
            throw std::invalid_argument("N must be > 0");
        }else{
            throw std::invalid_argument("N must be < MAXNP1");
        }
    }
    return m_r_to_N+((N-1)*m_nr_elements);
}
void RadInt::NewQuadrature(double eta, double P, unsigned int N, unsigned int lambda){
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
bool RadInt::CheckConvergence(){
    return (
            power(m_T[m_current_T]-m_T[m_current_T-1],2)
                            <=
            fabs(m_T[m_current_T]-m_T[m_current_T-2])*m_epsilon
            );
}
void RadInt::NextT(){
    long double T = 0.0;
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

RadInt::RadInt(double epsilon){
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
RadInt::~RadInt(){
    free(m_abscissas);
    free(m_weights);
    free(m_r_to_N);
}
long double RadInt::GetInt(double eta, double P, unsigned int N, unsigned int lambda){
    //NewQuadrature(eta, P, N, lambda);
    //for (int i=0; i<MINNRINTS; ++i){
    //    NextT();
    //}
    ////NextT();
    ////NextT();
    //while (not(CheckConvergence()) && m_current_T<NRSTEPS-1 && m_T[m_current_T]>RADINTTOLERANCE){
    //    NextT();
    //}
    //if (not(CheckConvergence())){
    //    if (m_T[m_current_T]>RADINTTOLERANCE){
    //        std::cerr << std::endl << eta << " " << P << " " << N << " " << lambda <<
    //        std::endl << m_current_T << " " << m_T[m_current_T] << " " << m_T[m_current_T-1] << " " << m_T[m_current_T-2] << std::endl <<
    //        power(m_T[m_current_T]-m_T[m_current_T-1],2)
    //                        << " <= " <<
    //        fabs(m_T[m_current_T]-m_T[m_current_T-2]) << " * " << m_epsilon << std::endl;
    //        throw std::runtime_error("Radial quadrature did not converge.");
    //    }else{
    //        m_T[m_current_T] = 0.0;
    //    }
    //}
    //return m_T[m_current_T];
//    1/4 Sqrt[\[Pi]] \[Eta]^(1/2 (-1-NN-\[Lambda])) (P \[Eta])^\[Lambda]
//        
//        Gamma[1/2 (1+NN+\[Lambda])]
//        
//        Hypergeometric1F1Regularized[1/2 (2-NN+\[Lambda]),3/2+\[Lambda],-P^2 \[Eta]]
//
//double gsl_sf_hyperg_1F1(double a, double b, double x);
//double gsl_sf_gammainv(const double x);
//double gsl_sf_gamma(const double x);

    double etaP = eta*P;
    double a,b,z;
    a=0.5*(2-(int)N+(int)lambda);
    b=1.5 /*3.0/2.0*/ + lambda;
    z=-etaP*P;
    //std::cout << N << " " << lambda << a << " " << b << " " << eta << " " << P << std::endl << std::flush;
    double hyper = gsl_sf_hyperg_1F1(a,b,z) * gsl_sf_gammainv(b);

    //std::cout << "gamma" << std::endl << std::flush;
    double gammaarg = 0.5*(1+N+lambda);
    double gamma = gsl_sf_gamma(gammaarg);

    double result = 0.25 * sqrtPi * pow(eta,-gammaarg) * power(etaP,lambda) * gamma * hyper;
    //double result = hyper;
    //std::cerr << a << " " << b << " " << z << " " << hyper << std::flush;

    return result;
}
