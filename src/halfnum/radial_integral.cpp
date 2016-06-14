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

//#include <halfnum/kfunction.h>
#include <halfnum/radial_integral.h>

//#if (MINNRINTS<3)
//#error "MINNRINTS must be at least 3."
//#endif

#include <gsl/gsl_sf_hyperg.h>
//double gsl_sf_hyperg_1F1(double a, double b, double x);
#include <gsl/gsl_sf_gamma.h>
//double gsl_sf_gammainv(const double x);
//double gsl_sf_gamma(const double x);
#include <tgmath.h>

#include <iostream>

//static const double Pi     = acos(-1.0L);
static const double sqrtPi = sqrt(acos(-1.0L));
//static const double logof2 = log(2.0);

template <typename T>
inline T power(T base, unsigned int exp){
    T result = 1.0;
    while (exp){
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}
//inline long double power(long double base, unsigned int exp){
//    double result = 1.0;
//    while (exp){
//        if (exp & 1)
//            result *= base;
//        exp >>= 1;
//        base *= base;
//    }
//    return result;
//}
//inline int power(int base, unsigned int exp){
//    double result = 1.0;
//    while (exp){
//        if (exp & 1)
//            result *= base;
//        exp >>= 1;
//        base *= base;
//    }
//    return result;
//}
//
//void RadInt::FillAbscissasWeightsDifferentialsRadii(){
//    //BEWARE: this computes the abscissas for the successive quadratures in one array!!!
//    //Therefore, consider that the weights for the first quadratures need to be adjusted to be used with
//    //the last one!!!
//    //BEWARE: since the original quadrature integrates over x in [-1,1] but this one integrates
//    //over r in 0,inf the differential is different (dr/dx = pow(2,r-1)/Log[2]). This differential is
//    //computed here as well. // NOTE: this is no longer the case, the computation of the differentials, that is
//    double xk, wk, arg;
//    int p = MINP;
//    int count = 0;
//    for (int k=1; k<=p; ++k, ++count){
//        arg = k*Pi/(1.0*p+1.0);
//        double sinarg = sin(arg);
//        xk = (p+1-2*k)/(1.0*p+1.0) +
//            (2.0/Pi)*(1.0+2.0*sinarg*sinarg/3.0)*cos(arg)*sinarg;
//        m_abscissas[count] = xk;//1.0-log(1.0-xk)/logof2;
//        wk = 16.0/(3.0*(p+1))*sinarg*sinarg*sinarg*sinarg;
//        m_weights[count] = wk;// / ((1.0-xk)*logof2);
//    }
//    for (int step=1; step<NRSTEPS; ++step){
//        p = 2*p+1;
//        for (int k=1; k<=p; k+=2, ++count){
//            arg = k*Pi/(1.0*p+1.0);
//            double sinarg = sin(arg);
//            xk = (p+1-2*k)/(1.0*p+1.0) +
//                (2.0/Pi)*(1.0+2.0*sinarg*sinarg/3.0)*cos(arg)*sinarg;
//            m_abscissas[count] = xk;//log(2.0/(1.0-xk))/logof2;
//            wk = 16.0/(3.0*(p+1))*sinarg*sinarg*sinarg*sinarg;
//            m_weights[count] = wk;// / ((1.0-xk)*logof2);
//        }
//    }
//    //double* r_to_N_it = m_r_to_N;
//    //for (unsigned int N=1; N<MAXNP1; ++N){
//    //    for (unsigned long int rit=0; rit<m_nr_elements; ++rit, ++r_to_N_it){
//    //        *r_to_N_it = power(m_abscissas[rit],N);
//    //    }
//    //}
//}
////double* RadInt::GetRStart(unsigned int N) const {
////    if (N==0 || N>=MAXNP1){
////        if (N==0){
////            throw std::invalid_argument("N must be > 0");
////        }else{
////            throw std::invalid_argument("N must be < MAXNP1");
////        }
////    }
////    return m_r_to_N+((N-1)*m_nr_elements);
////}
//void RadInt::NewQuadrature(double eta, double P, unsigned int N, unsigned int lambda){
//    m_absit = m_abscissas;
//    m_weiit = m_weights;
//    //m_radit = GetRStart(N);
//    m_N     = N;
//    m_current_T = -1;
//    m_P   = P;
//    m_eta = eta;
//    m_p   = MINP;
//    m_lambda = lambda;
//    double rsqrteta = 3.03485/sqrt(eta);
//    //m_start  = m_P - rsqrteta;
//    //m_end    = m_P + rsqrteta;
//    if (rsqrteta>m_P){
//        m_scale  = 0.5*(m_P + rsqrteta);
//        m_offset = m_scale;
//    }else{
//        m_offset = m_P;
//        m_scale  = rsqrteta;
//    }
//    for (int i=0; i<NRSTEPS; ++i){
//        m_T[i] = 0.0;
//    }
//}
//bool RadInt::CheckConvergence(){
//    return (
//            power(m_T[m_current_T]-m_T[m_current_T-1],2)
//                            <=
//            fabs(m_T[m_current_T]-m_T[m_current_T-2])*m_epsilon
//            );
//}
//void RadInt::NextT(){
//    long double T = 0.0;
//    if (m_current_T>=0){
//        T = 0.5*m_T[m_current_T];
//    }
//    for (unsigned int k=1; k<=m_p; k+=2, ++m_absit, ++m_weiit){
//        double r = m_scale * (*m_absit) + m_offset;
//        T += (*m_weiit) * power(r,m_N) * Kfunction(2.0*m_eta*m_P*r,m_lambda) * exp(-m_eta*(r - m_P)*(r - m_P));
//    }
//    m_p = 2*m_p+1;
//    ++m_current_T;
//    m_T[m_current_T] = m_scale * T;
//}
//
//RadInt::RadInt(double epsilon){
//    m_epsilon       = epsilon;
//    m_nr_elements   = power(2,NRSTEPS)-1;
//    m_abscissas     = (double*) malloc(m_nr_elements*sizeof(double));
//    m_weights       = (double*) malloc(m_nr_elements*sizeof(double));
//    memset(m_abscissas    ,0.0,m_nr_elements*sizeof(double));
//    memset(m_weights      ,0.0,m_nr_elements*sizeof(double));
//    //m_r_to_N = (double*) malloc((MAXNP1-1)*m_nr_elements*sizeof(double));
//    //memset(m_r_to_N       ,0.0,(MAXNP1-1)*m_nr_elements*sizeof(double));
//    FillAbscissasWeightsDifferentialsRadii();
//}
//RadInt::~RadInt(){
//    free(m_abscissas);
//    free(m_weights);
//    free(m_r_to_N);
//}
//long double GetRadInt(double eta, double P, unsigned int N, unsigned int lambda){
//    //NewQuadrature(eta, P, N, lambda);
//    //for (int i=0; i<MINNRINTS; ++i){
//    //    NextT();
//    //}
//    ////NextT();
//    ////NextT();
//    //while (not(CheckConvergence()) && m_current_T<NRSTEPS-1 && m_T[m_current_T]>RADINTTOLERANCE){
//    //    NextT();
//    //}
//    //if (not(CheckConvergence())){
//    //    if (m_T[m_current_T]>RADINTTOLERANCE){
//    //        std::cerr << std::endl << eta << " " << P << " " << N << " " << lambda <<
//    //        std::endl << m_current_T << " " << m_T[m_current_T] << " " << m_T[m_current_T-1] << " " << m_T[m_current_T-2] << std::endl <<
//    //        power(m_T[m_current_T]-m_T[m_current_T-1],2)
//    //                        << " <= " <<
//    //        fabs(m_T[m_current_T]-m_T[m_current_T-2]) << " * " << m_epsilon << std::endl;
//    //        throw std::runtime_error("Radial quadrature did not converge.");
//    //    }else{
//    //        m_T[m_current_T] = 0.0;
//    //    }
//    //}
//    //return m_T[m_current_T];
////    1/4 Sqrt[\[Pi]] \[Eta]^(1/2 (-1-NN-\[Lambda])) (P \[Eta])^\[Lambda]
////        
////        Gamma[1/2 (1+NN+\[Lambda])]
////        
////        Hypergeometric1F1Regularized[1/2 (2-NN+\[Lambda]),3/2+\[Lambda],-P^2 \[Eta]]
////
////double gsl_sf_hyperg_1F1(double a, double b, double x);
////double gsl_sf_gammainv(const double x);
////double gsl_sf_gamma(const double x);
//
//    double etaP = eta*P;
//    double z =-etaP*P;
//    double result;
//    if (N==5 && lambda==0 && -z > 500){
//        //treat special case that causes overflows in GSL
//        double summand1   = exp(z)*(5.0+2.0*P*P)/(8.0*eta*eta);
//        double summand2_1 = (3.0-4.0*z*(3.0-z)) / (16.0*etaP*pow(eta,2.5));
//        double summand2_2 = sqrtPi * erf(sqrt(-z));
//        result = summand1 + summand2_1*summand2_2;
//    }else{
//        double a,b;
//        a=0.5*(2-(int)N+(int)lambda);
//        b=1.5 /*3.0/2.0*/ + lambda;
//        //std::cout << "        " << N << " " << lambda << " " << a << " " << b << " " << z << " " << P << std::endl << std::flush;
//        //std::cout << "                 hyper " << std::flush;
//        double hyper = gsl_sf_hyperg_1F1(a,b,z);
//        //std::cout << hyper << std::flush;
//        //std::cout << " invgamma " << std::flush;
//        hyper *= gsl_sf_gammainv(b);
//        //std::cout << hyper << std::flush;
//
//        //std::cout << "gamma" << std::endl << std::flush;
//        double gammaarg = 0.5*(1+N+lambda);
//        //std::cout << " gamma " << std::flush;
//        double gamma = gsl_sf_gamma(gammaarg);
//        //std::cout << gamma << std::flush;
//        //std::cout << " result " << std::endl << std::flush;
//
//        result = 0.25 * sqrtPi * pow(eta,-gammaarg) * power(etaP,lambda) * gamma * hyper;
//        //double result = hyper;
//        //std::cerr << a << " " << b << " " << z << " " << hyper << std::flush;
//    }
void RadInt::Init(double eta_in, double P_in){
    eta      = eta_in;
    P        = P_in;
    //etaP     = eta*P;
    etaPP    = eta*P*P;
    //etaeta   = eta*eta;
    PP       = P*P;
    double sqrteta  = sqrt(eta);
    //eta3half = sqrteta*eta;
    //etaPetaP = etaP*etaP;
    erfetaP  = sqrtPi*erf(sqrteta*P);
    expetaPP = exp(-etaPP);
    _eta      = 1.0/eta;
    _sqrteta  = 1.0/sqrteta;
    _P        = 1.0/P;
    _etaeta   = _eta*_eta;
    _eta3half = _eta * _sqrteta;
    _etaP     = _eta*_P;
    _etaPP    = _etaP*_P;
    _etaPetaP = _etaP*_etaP;
}
double RadInt::GetRadInt(
            //double sqrteta  ,
            //double eta3half ,
            //double erfetaP  ,
            //double etaP     ,
            //double etaPP    ,
            //double expetaPP ,
            //double eta      ,
            //double P        ,
            int N           ,
            int lambda
        ){
    double result = 0.0;
    //double sqrteta  = sqrt(eta);
    //double eta3half = sqrteta * eta;
    //double erfetaP  = erf(sqrteta*P);
    //double etaP     = eta*P;
    //double etaPP    = etaP*P;
    //double expetaPP = exp(-etaPP);
    switch (N-1){
        case 0:
            switch (lambda){
                case 0:
                    result = 0.25*erfetaP*_P*_eta3half;
                    break;
                default:
                    throw std::runtime_error("Wrong input value for lambda.");
                    break;
            }
            break;
        case 1:
            switch (lambda){
                case 0:
                    result = 0.25*sqrtPi*_eta3half;
                    break;
                case 1:
                    //result = (2*expetaPP*P*sqrteta + sqrtPi*(2*etaPP - 1)*erfetaP) / (8 * etaP * etaP * sqrteta);
                    result = 0.25*expetaPP*_etaPP + 0.125*erfetaP*_etaPP*_sqrteta + 0.25*erfetaP*_eta3half;
                    break;
                default:
                    throw std::runtime_error("Wrong input value for lambda.");
                    break;
            }
            break;
        case 2:
            switch (lambda){
                case 0:
                    result = 0.25*expetaPP*_etaeta + (2*etaPP+1)*erfetaP*0.125*_etaP*_eta3half;
                    break;
                case 1:
                    result = 0.25*sqrtPi*P*_eta3half;
                    break;
                case 2:
                    //result = ((2*expetaPP*P*sqrteta*(2*etaPP-3)) + sqrtPi*(3+4*etaPP*(etaPP-1))*erfetaP)/(16*etaP*etaP*etaP*sqrteta);
                    result = - 0.375*expetaPP*_etaPetaP*_eta
                             + 0.25*expetaPP*_etaeta 
                             + 0.1875*erfetaP*_etaPetaP*_etaP*_sqrteta 
                             - 0.25*erfetaP*_etaP*_eta3half
                             + 0.25*P*erfetaP*_eta3half;
                    break;
                default:
                    throw std::runtime_error("Wrong input value for lambda.");
                    break;
            }
            break;
        case 3:
            switch (lambda){
                case 0:
                    result = 0.125*sqrtPi*(2*etaPP+3)*_eta3half*_eta;
                    break;
                case 1:
                    //result = ((2*expetaPP*P*sqrteta*(2*etaPP+1)) + sqrtPi*(-1+4*etaPP*(etaPP+1))*erfetaP)/(16*etaP*etaP*etaP*sqrteta);
                    result =   0.125*expetaPP*_etaP*_etaeta
                             + 0.25*expetaPP*P*_etaeta
                             - 0.0625*erfetaP*_etaPetaP*_eta3half
                             + 0.25*erfetaP*_eta*_eta3half
                             + 0.25*PP*erfetaP*_eta3half;
                    break;
                case 2:
                    result = 0.25*PP*sqrtPi*_eta3half;
                    break;
                case 3:
                    //result = expetaPP*(2*P*sqrteta*(15-8*etaPP+4*etaPP*etaPP) + sqrtPi*(-15+18*etaPP-12*etaPP*etaPP+8*etaPP*etaPP*etaPP)*erfetaP)/(32*etaP*etaP*etaP*etaP*sqrteta);
                    result =   0.9375*expetaPP*_etaPetaP*_etaP*_eta
                             - 0.5*expetaPP*_etaP*_etaeta
                             + 0.25*expetaPP*P*_etaeta
                             - 0.46875*erfetaP*_etaPetaP*_etaPetaP*_sqrteta
                             + 0.5625*erfetaP*_etaPetaP*_eta3half
                             - 0.375*erfetaP*_eta*_eta3half
                             + 0.25*PP*erfetaP*_eta3half;
                    break;
                default:
                    throw std::runtime_error("Wrong input value for lambda.");
                    break;
            }
            break;
        case 4:
            switch (lambda){
                case 0:
                    //result = (expetaPP*(5+2*etaPP))/(6*eta3half*eta3half) + sqrtPi*(3+4*etaPP*(3+etaPP))*erfetaP/(16*etaP*eta*eta3half);
                    result =   0.625*expetaPP*_etaeta*_eta
                             + 0.25*expetaPP*PP*_etaeta
                             + 0.1875*erfetaP*_etaP*_etaeta*_sqrteta
                             + 0.75*P*erfetaP*_eta*_eta3half
                             + 0.25*PP*P*erfetaP*_eta3half;
                    break;
                case 1:
                    result = 0.125*P*sqrtPi*(5+2*etaPP)*_eta*_eta3half;
                    break;
                case 2:
                    //result = ((2*expetaPP*P*sqrteta*(-3+4*etaPP*(etaPP+1))) + sqrtPi*(3+2*etaPP*(-3+6*etaPP+4*etaPP*etaPP))*erfetaP)/(32*etaP*etaP*etaP*eta3half);
                    result = - 0.1875*expetaPP*_etaPetaP*_etaeta
                             + 0.25*expetaPP*_etaeta*_eta
                             + 0.25*expetaPP*PP*_etaeta
                             + 0.09375*erfetaP*_etaPetaP*_etaP*_eta3half
                             - 0.1875*erfetaP*_etaP*_eta3half*_eta
                             + 0.375*P*erfetaP*_eta*_eta3half
                             + 0.25*PP*P*erfetaP*_eta3half;
                    break;
                case 3:
                    result = 0.25*PP*P*sqrtPi*_eta3half;
                    break;
                case 4:
                    //{
                    //double temp = sqrtPi*(105+8*etaPP*(-15+etaPP*(9+2*etaPP*(-2+etaPP))))*erfetaP/(2*sqrteta);
                    //result = (expetaPP*P*(-105+2*etaPP*(25+2*etaPP*(-5+2*etaPP))) + temp)/(32*etaP*etaP*etaP*etaP*etaP);
                    //}
                    result = - 3.28125*expetaPP*_etaPetaP*_etaPetaP*_eta
                             + 1.5625*expetaPP*_etaPetaP*_etaeta
                             - 0.625*expetaPP*_etaeta*_eta
                             + 0.25*expetaPP*PP*_etaeta
                             + 1.640625*erfetaP*_etaPetaP*_etaPetaP*_etaP*_sqrteta
                             - 1.875*erfetaP*_etaPetaP*_etaP*_eta3half
                             + 1.125*erfetaP*_etaP*_eta3half*_eta
                             - 0.5*P*erfetaP*_eta*_eta3half
                             + 0.25*PP*P*erfetaP*_eta3half;
                    break;
                default:
                    throw std::runtime_error("Wrong input value for lambda.");
                    break;
            }
            break;
        default:
            throw std::runtime_error("Wrong input value for N.");
    }
    
    //std::cout << result << std::endl;
    return result;
}
