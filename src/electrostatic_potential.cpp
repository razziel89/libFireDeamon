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
#include <electrostatic_potential.h>

void* _potentialThread(void* data){
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED,NULL);
    struct timespec req = {0/*req.tv_sec*/, 1L/*req.tv_nsec*/};
    //req.tv_sec = 0;
    //req.tv_nsec = 1L;
    GPSubData<double>* dat = (GPSubData<double>*) data;
    double *pnts = dat->GetData(0);
    double *ccos = dat->GetData(1);
    double *pots = dat->GetDataOutput();
    const int nr_ccos = dat->GetNr(1);
    const int nr_pots = dat->GetNrOutput();
    const int progress = 25;
    const bool progress_reports = dat->GetProgressReports();
    int* progress_bar = dat->GetProgressBar();
    pthread_mutex_t* mut = dat->GetMutex();

    double* p = pots;
    double* pc = pnts;
    if ( progress_reports ){
        for(int ip=0; ip<nr_pots;){
            for (int prog=0; prog<progress && ip<nr_pots; ++prog, ++ip, ++p){
                double xpc = *pc++;
                double ypc = *pc++;
                double zpc = *pc++;
                double sum = 0.0;
                double* c = ccos;
                for(int ic=0; ic<nr_ccos; ic+=4){
                    double norm = 0.0;
                    norm += (*c-xpc)*(*c-xpc); ++c;
                    norm += (*c-ypc)*(*c-ypc); ++c;
                    norm += (*c-zpc)*(*c-zpc); ++c;
                    sum += *c/(sqrt(norm));
                    ++c;
                }
                *p = sum;
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
        for(int ip=0; ip<nr_pots; ++ip, ++p){
            double xpc = *pc++;
            double ypc = *pc++;
            double zpc = *pc++;
            double sum = 0.0;
            double* c = ccos;
            for(int ic=0; ic<nr_ccos; ic+=4){
                double norm = 0.0;
                norm += (*c-xpc)*(*c-xpc); ++c;
                norm += (*c-ypc)*(*c-ypc); ++c;
                norm += (*c-zpc)*(*c-zpc); ++c;
                sum += *c/(sqrt(norm));
                ++c;
            }
            *p = sum;
            nanosleep(&req, (struct timespec *)NULL);
        }
    }
    pthread_exit(NULL);
}

void electrostatic_potential (bool progress_reports, int num_points, std::vector<double> points, std::vector<double> charges_coordinates, std::vector<double> *potential)
{
    //initialize everything
    PG globals; 
    init_parallel_generic(&progress_reports, &globals);

    //reserve data structures and fill them with input
    std::vector< std::vector<double> > input;
    input.reserve(2);
    const int split_col = 0;
    const int split_factor = 3;
    input.push_back(points);
    input.push_back(charges_coordinates);
    
    //fill class that holds data for each thread
    GPData<double> *data;
    try
    {
        data = new GPData<double>(progress_reports, globals.nr_threads, input, potential, &(globals.mutex), &(globals.progress_bar), split_col, split_factor);
    }
    catch( const std::invalid_argument& e ) {
        throw;
    }
    //perform computation
    do_parallel_generic<double>(_potentialThread, &globals, progress_reports, num_points, data);
    //transfer output data
    data->TransferOutput();
    //clean up
    //fprintf(stderr,"\nB Length: %d\n",potential->size());fflush(stderr);
    delete data;
    //fprintf(stderr,"\nE Length: %d\n",potential->size());fflush(stderr);
    finalize_parallel_generic(progress_reports, &globals);
    pg_global = NULL;
}
