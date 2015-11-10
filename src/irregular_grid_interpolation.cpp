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
#include <limits>
#include <stdexcept>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <parallel_generic.h>
#include <irregular_grid_interpolation.h>

void* _nearestInterpolationThread(void* data){
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED,NULL);
    struct timespec req = {0/*req.tv_sec*/, 1L/*req.tv_nsec*/};
    //req.tv_sec = 0;
    //req.tv_nsec = 1L;
    GPSubData<double>* dat = (GPSubData<double>*) data;
    double *inpnts = dat->GetData(0); //interpolation points
    double *pnts   = dat->GetData(1);
    double *vals   = dat->GetData(2);
    double *interp = dat->GetDataOutput();
    const int nr_pnts   = dat->GetNr(1);
    const int nr_interp = dat->GetNrOutput();
    const int progress = 250;
    const bool progress_reports = dat->GetProgressReports();
    int* progress_bar = dat->GetProgressBar();
    pthread_mutex_t* mut = dat->GetMutex();

    double* p  = interp;
    double* pc = inpnts;
    if ( progress_reports ){
        for(int ip=0; ip<nr_interp;){
            for (int prog=0; prog<progress && ip<nr_interp; ++prog, ++ip, ++p){
                double min_dist_squ = std::numeric_limits<double>::max();
                double min_dist_val = 0.0;
                double xpc = *pc++;
                double ypc = *pc++;
                double zpc = *pc++;
                double* c  = pnts;
                double* v  = vals;
                for(int ic=0; ic<nr_pnts; ic+=3){
                    double normsquared = 0.0;
                    normsquared += (*c-xpc)*(*c-xpc); ++c;
                    normsquared += (*c-ypc)*(*c-ypc); ++c;
                    normsquared += (*c-zpc)*(*c-zpc); ++c;
                    if ( normsquared < min_dist_squ ){
                        min_dist_squ = normsquared;
                        min_dist_val = *v;
                    }
                    ++v;
                }
                *p = min_dist_val;
            }
            //the nanosleep function serves as a possible point where Ctrl-C
            //can interrupt the programme
            nanosleep(&req, (struct timespec *)NULL);
            //the programme should not be cancelled while the mutex is locked
            pthread_setcancelstate(PTHREAD_CANCEL_DISABLE,NULL);
            pthread_mutex_lock(mut);
            *progress_bar += progress;
            pthread_mutex_unlock(mut);
            pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
        }
    }
    else {
        for(int ip=0; ip<nr_interp; ++ip, ++p){
            double min_dist_squ = std::numeric_limits<double>::max();
            double min_dist_val = 0.0;
            double xpc = *pc++;
            double ypc = *pc++;
            double zpc = *pc++;
            double* c  = pnts;
            double* v  = vals;
            for(int ic=0; ic<nr_pnts; ic+=3){
                double normsquared = 0.0;
                normsquared += (*c-xpc)*(*c-xpc); ++c;
                normsquared += (*c-ypc)*(*c-ypc); ++c;
                normsquared += (*c-zpc)*(*c-zpc); ++c;
                if ( normsquared < min_dist_squ ){
                    min_dist_squ = normsquared;
                    min_dist_val = *v;
                }
                ++v;
            }
            *p = min_dist_val;
            nanosleep(&req, (struct timespec *)NULL);
        }
    }
    pthread_exit(NULL);
}

void* _inverseDistanceWeightingInterpolationThread(void* data){
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED,NULL);
    struct timespec req = {0/*req.tv_sec*/, 1L/*req.tv_nsec*/};
    //req.tv_sec = 0;
    //req.tv_nsec = 1L;
    GPSubData<double>* dat = (GPSubData<double>*) data;
    double *inpnts = dat->GetData(0); //interpolation points
    double *pnts   = dat->GetData(1);
    double *vals   = dat->GetData(2);
    double *config = dat->GetData(3);
    const double distance_exponent =      config[0];
    const int    distance_function = (int)config[1];
    const double combined_exponent = -distance_exponent*1.0/distance_function;
    double *interp = dat->GetDataOutput();
    const int nr_pnts   = dat->GetNr(1);
    const int nr_interp = dat->GetNrOutput();
    const int progress = 25;
    const bool progress_reports = dat->GetProgressReports();
    int* progress_bar = dat->GetProgressBar();
    pthread_mutex_t* mut = dat->GetMutex();
    
    double* p  = interp;
    double* pc = inpnts;
    if ( progress_reports ){
        for(int ip=0; ip<nr_interp;){
            for (int prog=0; prog<progress && ip<nr_interp; ++prog, ++ip, ++p){
                long double sum_weights = 0.0L;
                long double sum_weiths_times_value = 0.0L;
                double xpc = *pc++;
                double ypc = *pc++;
                double zpc = *pc++;
                double* c  = pnts;
                double* v  = vals;
                for(int ic=0; ic<nr_pnts; ic+=3){
                    double weight = 0.0;
                    double t, temp;
                    t = fabs(*c-xpc); ++c;
                    temp = t;
                    for (int j=1; j<distance_function; ++j){
                        temp *= t;
                    }
                    weight += temp;
                    t = fabs(*c-ypc); ++c;
                    temp = t;
                    for (int j=1; j<distance_function; ++j){
                        temp *= t;
                    }
                    weight += temp;
                    t = fabs(*c-zpc); ++c;
                    temp = t;
                    for (int j=1; j<distance_function; ++j){
                        temp *= t;
                    }
                    weight += temp;
                    //weight = pow(weight,1.0/distance_function);
                    //weight = 1.0/pow(weight,distance_exponent);
                    //weight = pow(weight,-distance_exponent*1.0/distance_function);
                    weight = pow(weight,combined_exponent);
                    sum_weights += (long double)weight;
                    sum_weiths_times_value += (long double)(weight * (*v));
                    ++v;
                }
                *p = (double)(sum_weiths_times_value/sum_weights);
//                fprintf(stdout,"Exp: %f, Func: %d, Weight: %Lf, Weight*Val: %Lf, Quot: %Lf, DiffL %f\n",distance_exponent,distance_function,sum_weights,sum_weiths_times_value,sum_weiths_times_value/sum_weights,((double)(sum_weiths_times_value/sum_weights))-min_dist_val);
            }
//            fflush(stdout);
            //the nanosleep function serves as a possible point where Ctrl-C
            //can interrupt the programme
            nanosleep(&req, (struct timespec *)NULL);
            //the programme should not be cancelled while the mutex is locked
            pthread_setcancelstate(PTHREAD_CANCEL_DISABLE,NULL);
            pthread_mutex_lock(mut);
            *progress_bar += progress;
            pthread_mutex_unlock(mut);
            pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
        }
    }
    else {
        for(int ip=0; ip<nr_interp; ++ip, ++p){
            long double sum_weights = 0.0L;
            long double sum_weiths_times_value = 0.0L;
            double xpc = *pc++;
            double ypc = *pc++;
            double zpc = *pc++;
            double* c  = pnts;
            double* v  = vals;
            for(int ic=0; ic<nr_pnts; ic+=3){
                double weight = 0.0;
                double t, temp;
                t = fabs(*c-xpc); ++c;
                temp = t;
                for (int j=1; j<distance_function; ++j){
                    temp *= t;
                }
                weight += temp;
                t = fabs(*c-ypc); ++c;
                temp = t;
                for (int j=1; j<distance_function; ++j){
                    temp *= t;
                }
                weight += temp;
                t = fabs(*c-zpc); ++c;
                temp = t;
                for (int j=1; j<distance_function; ++j){
                    temp *= t;
                }
                weight += temp;
                weight = pow(weight,distance_exponent*1.0/distance_function);
                sum_weights += (long double)weight;
                sum_weiths_times_value += (long double)(weight* (*v));
                ++v;
            }
            *p = (double)(sum_weiths_times_value/sum_weights);
            nanosleep(&req, (struct timespec *)NULL);
        }
    }
    pthread_exit(NULL);
}

void generic_interpolation(bool progress_reports, int num_interpolation_points, std::vector<double> points, std::vector<double> values, std::vector<double> interpolation_points, std::vector<double> *interpolation, int interpolation_type, double distance_exponent, int distance_function)
{
    void* (*interp_func)(void*);
    //initialize everything
    PG globals; 
    init_parallel_generic(&progress_reports, &globals);

    //reserve data structures and fill them with input
    std::vector< std::vector<double> > input;
    input.reserve(3);
    const int split_col = 0;
    const int split_factor = 3;
    input.push_back(interpolation_points);
    input.push_back(points);
    input.push_back(values);
    
    //fill class that holds data for each thread
    GPData<double> *data;
    std::vector<double> config;
    switch (interpolation_type){
        case 1:
            if (progress_reports){
                fprintf(stdout,"Using nearest-neighbour interpolation.\n");
            }
            interp_func = _nearestInterpolationThread;
            try
            {
                data = new GPData<double>(progress_reports, globals.nr_threads, input, interpolation, &(globals.mutex), &(globals.progress_bar), split_col, split_factor);
            }
            catch( const std::invalid_argument& e ) {
                throw;
            }
            break;
        case 2:
            assert( distance_exponent >= 0.0 && distance_function >= 1);
            if (progress_reports){
                fprintf(stdout,"Using inverse distance weighting with p=%f and the ",distance_exponent);
                switch (distance_function%10){
                    case 1:
                        fprintf(stdout,"%dst",distance_function);
                        break;
                    case 2:
                        fprintf(stdout,"%dnd",distance_function);
                        break;
                    case 3:
                        fprintf(stdout,"%drd",distance_function);
                        break;
                    default:
                        fprintf(stdout,"%dth",distance_function);
                        break;
                }
                fprintf(stdout," root.\n");
            }
            input.reserve(input.size()+1);
            config.reserve(2);
            config.push_back(distance_exponent);
            config.push_back(distance_function);
            input.push_back(config);
            interp_func = _inverseDistanceWeightingInterpolationThread;
            try
            {
                data = new GPData<double>(progress_reports, globals.nr_threads, input, interpolation, &(globals.mutex), &(globals.progress_bar), split_col, split_factor);
            }
            catch( const std::invalid_argument& e ) {
                throw;
            }
            break;
        default:
            throw std::invalid_argument( "Supported interpolationmethods: 1: nearest neighbour. Choose one of them." );
            break;
    }
    fflush(stdout);
    //perform computation
    do_parallel_generic<double>(interp_func, &globals, progress_reports, num_interpolation_points, data);
    //transfer output data
    data->TransferOutput();
    //clean up
    delete data;
    finalize_parallel_generic(progress_reports, &globals);
    pg_global = NULL;
}
