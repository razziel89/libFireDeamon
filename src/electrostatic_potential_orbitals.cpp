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
#include <cstdlib>
#include <pthread.h>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <parallel_generic.h>
#include <halfnum/radial_integral.h>
#include <halfnum/angular_integral.h>
#include <electrostatic_potential_orbitals.h>
#include <gsl/gsl_sf_legendre.h>

#include <iostream>

static const double Pi    = acos(-1.0L);
static const double sqrt2 = sqrt(2.0);

inline double spherical_harmonic(int l, int m, double theta, double phi){
    int abs_m = m ? m>=0 : -m;
    double result = gsl_sf_legendre_sphPlm (l, abs_m, cos(theta));
    if     (m > 0){
        result *= sqrt2 * (((abs_m+1) % 2)-((abs_m) % 2)) * cos(abs_m*phi);
    }
    else{if(m < 0){
        result *= sqrt2 * (((abs_m+1) % 2)-((abs_m) % 2)) * sin(abs_m*phi);
    }}
    return result;
}

int factorial[] = {
    1                  ,
    1                  ,
    2                  ,
    6                  ,
    24                 ,
    120                ,
    720                ,
    5040               ,
    40320              ,
    362880             ,
    3628800            ,
    39916800           ,
    479001600
};

//int dfactorial[] = {
//    1                               ,
//    1                               ,
//    2                               ,
//    3                               ,
//    8                               ,
//    15                              ,
//    48                              ,
//    105                             ,
//    384                             ,
//    945                             ,
//    3840                            ,
//    10395                           ,
//    46080                           ,
//    135135                          ,
//    645120                          ,
//    2027025                         ,
//    10321920                        ,
//    34459425                        ,
//    185794560                       ,
//    654729075                       ,
//    3715891200//                      ,
//}

inline int binomial(int n, int k){
    return (factorial[n])/(factorial[n-k]*factorial[k]);
}

inline void cartesian_to_spherical (double x, double y, double z, double &r, double &theta, double &phi){
    r = sqrt(x*x + y*y + z*z);
    theta=acos(z/r);
    phi=atan2(y, x);
    return;
}

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

inline double X_mu_nu(double A[3], double B[3], double P[3], unsigned int a[3], unsigned int b[3], double d_mu_nu, double eta_mu_nu, RadInt& RI, const AngInt& AI){
    double result = 0.0;
    double absP;
    double thetaP;
    double phiP;
    cartesian_to_spherical(P[0],P[1],P[2],absP,thetaP,phiP);
    //loop over all ai
    //int a_minus1 = ((a[0]+a[1]+a[2])%2==0) ? 1 : -1;
    for (unsigned int ax=0; ax<=a[0]; ++ax/*, a_minus1=-a_minus1*/){
        int a_binom[3];
        a_binom[0] = binomial(a[0],ax);
    for (unsigned int ay=0; ay<=a[1]; ++ay/*, a_minus1=-a_minus1*/){
        a_binom[1] = binomial(a[1],ay) * a_binom[0];
    for (unsigned int az=0; az<=a[2]; ++az/*, a_minus1=-a_minus1*/){
        a_binom[2] = binomial(a[2],az) * a_binom[1];

        int a_minus1 = ((a[0]+a[1]+a[2]+ax+ay+az)%2==0) ? 1 : -1;
        unsigned int sum_a = ax+ay+az;
        double Apow = d_mu_nu * a_binom[2] * a_minus1 * power(A[0],a[0]-ax) * power(A[1],a[1]-ay) * power(A[2],a[2]-az);

        //loop over all bi
        //int b_minus1 = ((b[0]+b[1]+b[2])%2==0) ? 1 : -1;
        for (unsigned int bx=0; bx<=b[0]; ++bx/*, b_minus1=-b_minus1*/){
            int b_binom[3];
            b_binom[0] = binomial(b[0],bx);
        for (unsigned int by=0; by<=b[1]; ++by/*, b_minus1=-b_minus1*/){
            b_binom[1] = binomial(b[1],by) * b_binom[0];
        for (unsigned int bz=0; bz<=b[2]; ++bz/*, b_minus1=-b_minus1*/){
            b_binom[2] = binomial(b[2],bz) * b_binom[1];

            int b_minus1 = ((b[0]+b[1]+b[2]+bx+by+bz)%2==0) ? 1 : -1;
            unsigned int sum_b = bx+by+bz;
            double Bpow = b_binom[2] * b_minus1 * power(B[0],b[0]-bx) * power(B[1],b[1]-by) * power(B[2],b[2]-bz);
            
            double Q = 0.0;
            for (unsigned int lambda=0; lambda<=sum_a+sum_b; ++lambda){
                double Omega = 0.0;
                for (int gamma=-((int)lambda); gamma<=(int)(lambda); ++gamma){
                    double integral = AI.GetInt(lambda, gamma, ax+bx, ay+by, az+bz);
                    if (integral>0.0){
                        Omega += spherical_harmonic(lambda, gamma, thetaP, phiP) * integral;
                    }
                }
                if (Omega>0.0){
                    Q += RI.GetInt(eta_mu_nu, absP, sum_a+sum_b+1, lambda) * Omega;
                }
            }
            result += Apow * Bpow * Q; 
        }}}
    }}}
    return result * 2.0*2.0*Pi;
}

void* _potentialThreadOrbitals(void* data){
    //std::cout << "thread0" << std::endl << std::flush;

    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED,NULL);
    struct timespec req = {0/*req.tv_sec*/, 1L/*req.tv_nsec*/};
    //req.tv_sec = 0;
    //req.tv_nsec = 1L;
    GPSubData<double,double,double,int,double,double,double,int,double>* dat;
    dat = static_cast<GPSubData<double,double,double,int,double,double,double,int,double>*>(data);
    //std::cout << "thread1" << std::endl << std::flush;

    //input = std::make_tuple(potential_grid,prim_centers,prim_angular,prim_exponents,prim_coefficients,P_matrix);
    double *grdpnts  = dat->GetData<0>(); //gridpoints
    double *prm_cent = dat->GetData<1>(); //centers of primitives
    int    *prm_ang  = dat->GetData<2>(); //angular exponents of primitives
    double *prm_exp  = dat->GetData<3>(); //exponents of primitives
    double *prm_coef = dat->GetData<4>(); //contraction coefficients of primitives
    double *P_mu_nu  = dat->GetData<5>(); //first order density matrix P_mu_nu
    int    *screen   = dat->GetData<6>(); //whether or not the overlap integral is sufficient
    double *pot      = dat->GetDataOutput();
    double *cutoff   = dat->GetData<7>();
    double cutoff_2  = (*cutoff)*(*cutoff);
    //std::cout << "thread2" << std::endl << std::flush;

    const int nr_prim  = dat->GetNr<1>()/3; //number of primitive gaussian basis functions
    const int nr_pot   = dat->GetNrOutput();
    const int progress = 25;
    const bool progress_reports = dat->GetProgressReports();
    int* progress_bar = dat->GetProgressBar();
    pthread_mutex_t* mut = dat->GetMutex();
    //std::cout << "thread3" << std::endl << std::flush;

    RadInt RI(1.0e-6);
    const AngInt AI;

    //std::cout << "thread init" << std::endl << std::flush;

    double* gp        = grdpnts;
    double* potential = pot;
    for(int pp = 0; pp < nr_pot;){
        for (int prog=0; prog<progress && pp<nr_pot; ++prog, ++pp){
            //std::cout << "point " << pp << "/" << nr_pot << std::endl << std::flush;
            double sum = 0.0;
            double r[3];
            r[0] = *gp; ++gp;
            r[1] = *gp; ++gp;
            r[2] = *gp; ++gp;
            //mu
            double* pcent_mu = prm_cent;
            int*    pang_mu  = prm_ang;
            double* pexp_mu  = prm_exp;
            double* pcoef_mu = prm_coef;
            //std::cout << "thread inner start" << std::endl << std::flush;
            for (int mu=0; mu<nr_prim; ++mu){
                double A[3];
                A[0] = *pcent_mu - r[0]; ++pcent_mu;
                A[1] = *pcent_mu - r[1]; ++pcent_mu;
                A[2] = *pcent_mu - r[2]; ++pcent_mu;
                unsigned int a[3];
                a[0] = (unsigned int)*pang_mu; ++pang_mu;
                a[1] = (unsigned int)*pang_mu; ++pang_mu;
                a[2] = (unsigned int)*pang_mu; ++pang_mu;
                double xi_mu = *pexp_mu; ++pexp_mu;
                double d_mu  = *pcoef_mu; ++pcoef_mu;
                //screening
                if ( ((A[0]*A[0]) + (A[1]*A[1]) + (A[2]*A[2])) < cutoff_2){

                    sum += P_mu_nu[mu*nr_prim + mu] * X_mu_nu(A,A,A,a,a,d_mu*d_mu,2.0*xi_mu,RI,AI);
                    //nu
                    //double* pcent_nu = prm_cent;
                    //int*    pang_nu  = prm_ang;
                    //double* pexp_nu  = prm_exp;
                    //double* pcoef_nu = prm_coef;
                    double* pcent_nu = pcent_mu;
                    int*    pang_nu  = pang_mu;
                    double* pexp_nu  = pexp_mu;
                    double* pcoef_nu = pcoef_mu;
                    for (int nu=mu+1; nu<nr_prim; ++nu){ //both, P_mu_nu and X_mu_nu are symmetric
                        double B[3];
                        B[0] = *pcent_nu - r[0]; ++pcent_nu;
                        B[1] = *pcent_nu - r[1]; ++pcent_nu;
                        B[2] = *pcent_nu - r[2]; ++pcent_nu;
                        unsigned int b[3];
                        b[0] = (unsigned int)*pang_nu; ++pang_nu;
                        b[1] = (unsigned int)*pang_nu; ++pang_nu;
                        b[2] = (unsigned int)*pang_nu; ++pang_nu;
                        double xi_nu = *pexp_nu; ++pexp_nu;
                        double d_nu  = *pcoef_nu; ++pcoef_nu;
                        //screening
                        if ( (B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) < cutoff_2){
                            //combined values
                            //if ( (P[0]*P[0] + P[1]*P[1] + P[2]*P[2]) < cutoff_2){
                            if ( screen[mu*nr_prim + nu] != 0 ){
                                double eta_mu_nu = xi_mu + xi_nu;
                                double P[3];
                                P[0] = (A[0]-B[0]);
                                P[1] = (A[1]-B[1]);
                                P[2] = (A[2]-B[2]);
                                double C_mu_nu = exp(-(d_mu*d_nu*(
                                                                    (P[0]*P[0])+(P[1]*P[1])+(P[2]*P[2])
                                                                 ))/eta_mu_nu);
                                P[0] = (xi_mu*A[0] + xi_nu*B[0]) / eta_mu_nu;
                                P[1] = (xi_mu*A[1] + xi_nu*B[1]) / eta_mu_nu;
                                P[2] = (xi_mu*A[2] + xi_nu*B[2]) / eta_mu_nu;
                                double d_mu_nu = d_mu * d_nu * C_mu_nu;
                                sum += 2.0 * P_mu_nu[mu*nr_prim + nu] * X_mu_nu(A,B,P,a,b,d_mu_nu,eta_mu_nu,RI,AI);
                            }
                        }
                    }
                }
            }
            //std::cout << "thread inner asign" << std::endl << std::flush;
            *potential = sum; ++potential;
            //std::cout << "thread inner end" << std::endl << std::flush;
        }
        //std::cout << "thread done" << std::endl << std::flush;
        //the nanosleep function serves as a possible point where Ctrl-C
        //can interrupt the programme
        nanosleep(&req, (struct timespec *)NULL);
        //the programme should not be cancelled while the mutex is locked
        if (progress_reports){
            pthread_setcancelstate(PTHREAD_CANCEL_DISABLE,NULL);
            pthread_mutex_lock(mut);
            *progress_bar += progress;
            pthread_mutex_unlock(mut);
            pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
        }
    }
    pthread_exit(NULL);
}

void electrostatic_potential_orbitals(bool progress_reports, int num_gridpoints, int num_primitives, std::vector<double> prim_centers, std::vector<double> prim_exponents, std::vector<double> prim_coefficients, std::vector<int> prim_angular, std::vector<double> potential_grid, std::vector<double> P_matrix, std::vector<int> screen, std::vector<double> *potential, double cutoff)
{
    //std::cout << "start" << std::endl << std::flush;
    //std::cout << prim_centers.size() << " " << prim_exponents.size() << " " << prim_coefficients.size() << " " << prim_angular.size() << " " << potential_grid.size() << " " << P_matrix.size() << " " << std::endl << std::flush;
    //initialize everything
    PG globals; 
    init_parallel_generic(&progress_reports, &globals);

    const int split_factor = 3;

    std::vector<double> cutoff_vec;
    cutoff_vec.push_back(cutoff);

    //reserve data structures and fill them with input
    tuple_of_vectors<double,double,int,double,double,double,int,double> input;
    input = std::make_tuple(potential_grid,prim_centers,prim_angular,prim_exponents,prim_coefficients,P_matrix,screen,cutoff_vec);
    //std::cout << "tuple" << std::endl << std::flush;

    //fill class that holds data for each thread
    GPData<double,double,double,int,double,double,double,int,double> *data;
    try
    {
        data = new GPData<double,double,double,int,double,double,double,int,double>(progress_reports, globals.nr_threads, input, potential, &(globals.mutex), &(globals.progress_bar), split_factor, 1, /*interlace=*/true);
    }
    catch( const std::invalid_argument& e ) {
        throw;
    }
    fflush(stdout);
    //perform computation
    //std::cout << "start do" << std::endl << std::flush;
    do_parallel_generic<double,double,double,int,double,double,double,int,double>(_potentialThreadOrbitals, &globals, progress_reports, num_gridpoints, data);
    //std::cout << "end do" << std::endl << std::flush;
    //transfer output data
    data->TransferOutput();
    //clean up
    delete data;
    finalize_parallel_generic(progress_reports, &globals);
    pg_global = NULL;
    //std::cout << "end" << std::endl << std::flush;
}
