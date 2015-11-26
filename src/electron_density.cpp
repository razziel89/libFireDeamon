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
#include <stdexcept>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <parallel_generic.h>
#include <electron_density.h>

void* _electronDensityThreadCutoff(void* data){
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED,NULL);
    struct timespec req = {0/*req.tv_sec*/, 1L/*req.tv_nsec*/};
    //req.tv_sec = 0;
    //req.tv_nsec = 1L;
    GPSubData<double,double>* dat = (GPSubData<double,double>*) data;

    double *grdpnts  = dat->GetData(0); //gridpoints
    double *prm_cent = dat->GetData(1); //centers of primitives
    double *prm_ang  = dat->GetData(2); //angular exponents of primitives
    double *prm_exp  = dat->GetData(3); //exponents of primitives
    double *prm_coef = dat->GetData(4); //contraction coefficients of primitives
    double *mo_coef  = dat->GetData(5); //molecular orbital coefficients
    double *cutoff   = dat->GetData(6); //the cutoff above which no density will be computed
    double cutoff_2  = (*cutoff)*(*cutoff); //only the square is needed
    double *dens     = dat->GetDataOutput();

    //const int nr_grdpnts = dat->GetNr(0);
    //const int nr_cent    = dat->GetNr(1);
    //const int nr_ang     = dat->GetNr(2);
    //const int nr_exp     = dat->GetNr(3);
    const int nr_coef    = dat->GetNr(4);
    const int nr_dens    = dat->GetNrOutput();
    //const int sub_nr = dat->GetSubNr();
    const int progress = 5000;
    const bool progress_reports = dat->GetProgressReports();
    int* progress_bar = dat->GetProgressBar();
    pthread_mutex_t* mut = dat->GetMutex();

    double* gp      = grdpnts;
    double* density = dens;
    for(int dp = 0; dp < nr_dens;){
        for (int prog=0; prog<progress && dp<nr_dens; ++prog, ++dp){
            double sum = 0.0;
            double xgp = *gp; ++gp;
            double ygp = *gp; ++gp;
            double zgp = *gp; ++gp;
            double* pcent = prm_cent;
            double* pang  = prm_ang;
            double* pexp  = prm_exp;
            double* pcoef = prm_coef;
            double* mo    = mo_coef;
            for(int ip =0; ip<nr_coef; ++ip){
                double dx = xgp-(*(pcent+0));
                double dy = ygp-(*(pcent+1));
                double dz = zgp-(*(pcent+2));
                double dist_2 = (dx*dx + dy*dy + dz*dz);
                if (dist_2 < cutoff_2){
                    double exp_factor = exp( -(*pexp) * dist_2 );
                    double prefac;
                    prefac  = pow(dx,(int)(*(pang+0)));
                    prefac *= pow(dy,(int)(*(pang+1)));
                    prefac *= pow(dz,(int)(*(pang+2)));
                    sum += (*pcoef) * (*mo) * prefac * exp_factor;
                }
                ++pcoef; ++mo; ++pexp; pang+=3; pcent+=3;
            }
            *density = sum*sum; ++density;
        }
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

void* _electronDensityThreadNoCutoff(void* data){
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED,NULL);
    struct timespec req = {0/*req.tv_sec*/, 1L/*req.tv_nsec*/};
    //req.tv_sec = 0;
    //req.tv_nsec = 1L;
    GPSubData<double,double>* dat = (GPSubData<double,double>*) data;

    double *grdpnts  = dat->GetData(0); //gridpoints
    double *prm_cent = dat->GetData(1); //centers of primitives
    double *prm_ang  = dat->GetData(2); //angular exponents of primitives
    double *prm_exp  = dat->GetData(3); //exponents of primitives
    double *prm_coef = dat->GetData(4); //contraction coefficients of primitives
    double *mo_coef  = dat->GetData(5); //molecular orbital coefficients
    double *dens     = dat->GetDataOutput();

    //const int nr_grdpnts = dat->GetNr(0);
    //const int nr_cent    = dat->GetNr(1);
    //const int nr_ang     = dat->GetNr(2);
    //const int nr_exp     = dat->GetNr(3);
    const int nr_coef    = dat->GetNr(4);
    const int nr_dens    = dat->GetNrOutput();
    //const int sub_nr = dat->GetSubNr();
    const int progress = 5000;
    const bool progress_reports = dat->GetProgressReports();
    int* progress_bar = dat->GetProgressBar();
    pthread_mutex_t* mut = dat->GetMutex();

    double* gp      = grdpnts;
    double* density = dens;
    for(int dp = 0; dp < nr_dens;){
        for (int prog=0; prog<progress && dp<nr_dens; ++prog, ++dp){
            double sum = 0.0;
            double xgp = *gp; ++gp;
            double ygp = *gp; ++gp;
            double zgp = *gp; ++gp;
            double* pcent = prm_cent;
            double* pang  = prm_ang;
            double* pexp  = prm_exp;
            double* pcoef = prm_coef;
            double* mo    = mo_coef;
            for(int ip =0; ip<nr_coef; ++ip){
                double dx = xgp-(*pcent); ++pcent;
                double dy = ygp-(*pcent); ++pcent;
                double dz = zgp-(*pcent); ++pcent;
                double dist_2 = (dx*dx + dy*dy + dz*dz);
                double exp_factor = exp( -(*pexp) * dist_2 );
                double prefac;
                prefac  = pow(dx,(int)(*(pang+0)));
                prefac *= pow(dy,(int)(*(pang+1)));
                prefac *= pow(dz,(int)(*(pang+2)));
                sum += (*pcoef) * (*mo) * prefac * exp_factor;
                ++pcoef; ++mo; ++pexp; pang+=3;
            }
            *density = sum*sum; ++density;
        }
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

void electron_density(bool progress_reports, int num_gridpoints, std::vector<double> prim_centers, std::vector<double> prim_exponents, std::vector<double> prim_coefficients, std::vector<double> prim_angular, std::vector<double> density_grid, std::vector<double> mo_coefficients, std::vector<double> *density, double cutoff)
{
    //initialize everything
    PG globals; 
    init_parallel_generic(&progress_reports, &globals);

    //check whether a cutoff shall be considered
    bool use_cutoff = (cutoff > 0.0);

    //reserve data structures and fill them with input
    std::vector< std::vector<double> > input;
    input.reserve(6+(use_cutoff ? 1 : 0));
    const int split_col = 0;
    const int split_factor = 3;
    input.push_back(density_grid);
    input.push_back(prim_centers);
    input.push_back(prim_angular);
    input.push_back(prim_exponents);
    input.push_back(prim_coefficients);
    input.push_back(mo_coefficients);

    std::vector<double> vec_cutoff;
    if (use_cutoff){
        vec_cutoff.push_back(cutoff);
        input.push_back(vec_cutoff);
    }
    
    //fill class that holds data for each thread
    GPData<double,double> *data;
    try
    {
        data = new GPData<double,double>(progress_reports, globals.nr_threads, input, density, &(globals.mutex), &(globals.progress_bar), split_col, split_factor, 1, /*interlace=*/true);
    }
    catch( const std::invalid_argument& e ) {
        throw;
    }
    fflush(stdout);
    //perform computation
    if (use_cutoff){
        do_parallel_generic<double,double>(_electronDensityThreadCutoff, &globals, progress_reports, num_gridpoints, data);
    }
    else{
        do_parallel_generic<double,double>(_electronDensityThreadNoCutoff, &globals, progress_reports, num_gridpoints, data);
    }
    //transfer output data
    data->TransferOutput();
    //clean up
    delete data;
    finalize_parallel_generic(progress_reports, &globals);
    pg_global = NULL;
}
