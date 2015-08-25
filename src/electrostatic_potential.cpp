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
#include <signal.h>
#include <time.h>

pthread_t* threads;
pthread_mutex_t mutex;
int nr_threads = 0;

void signal_callback_handler(int signum)
{
    fprintf(stderr, "\nCaught signal %d\n", signum);
    for ( int thread_count=0; thread_count<nr_threads; ++thread_count ){
        int rc = pthread_cancel(*(threads+thread_count));
        if(rc){
            fprintf(stderr, "Cancellation of thread number %d out of %d could not be initialized.\n", thread_count+1, nr_threads);
        }
    }
    delete[] threads;
    exit(signum);
}

//BEGINNING of class definitions
//this class holds all the data that will be worked on
//I have chosen to use a class in order to be able to pass a single argument to the thread functions
class SubPotData {

    public:

        SubPotData();
        SubPotData(bool progress_reports, int sub_nr, int l_pnts, int l_ccos, int l_pots, double* points, double* charges_coordinates, double* potential, pthread_mutex_t* mutex, int* progress_bar);
        ~SubPotData();
        double* GetPnts();
        double* GetCCos();
        double* GetPots();
        int GetNrPnts();
        int GetNrCCos();
        int GetNrPots();
        int GetSubNr();
        int* GetProgressBar();
        bool GetProgressReports();
        pthread_mutex_t* GetMutex();

    protected:

        double *m_pnts;
        double *m_ccos;
        double *m_pots;
        int m_nr_pnts;
        int m_nr_ccos;
        int m_nr_pots;
        int m_sub_nr;
        int* m_progress_bar;
        pthread_mutex_t* m_mut;
        bool m_progress_reports;

};

SubPotData::SubPotData(){
    m_pnts            = NULL;
    m_ccos            = NULL;
    m_pots            = NULL;
    m_nr_pnts         = 0;
    m_nr_ccos         = 0;
    m_nr_pots         = 0;
    m_sub_nr          = 0;
    m_progress_bar    = NULL;
    m_mut             = NULL;
    m_progress_reports= false;
}

SubPotData::SubPotData(bool progress_reports, int sub_nr, int l_pnts, int l_ccos, int l_pots, double* points, double* charges_coordinates, double* potential, pthread_mutex_t* mutex, int* progress_bar){
    m_progress_reports = progress_reports;
    m_mut              = mutex;
    m_nr_pnts = l_pnts;
    m_nr_ccos = l_ccos;
    m_nr_pots = l_pots;
    m_sub_nr = sub_nr;
    m_progress_bar = progress_bar;
    m_pnts = points;
    m_ccos = charges_coordinates;
    m_pots = potential;
}

SubPotData::~SubPotData(){}

double* SubPotData::GetPnts(){
    return m_pnts;
}
double* SubPotData::GetCCos(){
    return m_ccos;
}
double* SubPotData::GetPots(){
    return m_pots;
}
int SubPotData::GetNrPnts(){
    return m_nr_pnts;
}
int SubPotData::GetNrCCos(){
    return m_nr_ccos;
}
int SubPotData::GetNrPots(){
    return m_nr_pots;
}
int SubPotData::GetSubNr(){
    return m_sub_nr;
}
int* SubPotData::GetProgressBar(){
    return m_progress_bar;
}
bool SubPotData::GetProgressReports(){
    return m_progress_reports;
}
pthread_mutex_t* SubPotData::GetMutex(){
    return m_mut;
}

class PotData : public SubPotData {

    public:

        PotData();
        PotData(bool progress_reports, int nr_subs, std::vector<double> *points, std::vector<double> *charges_coordinates, std::vector<double> *potential, pthread_mutex_t* mutex, int* progress_bar);
        ~PotData();
        SubPotData* GetSubPotData(int index);
        void TransferPotential(std::vector<double> *pot);

    private:
        
        int m_nr_subs;
        std::vector<SubPotData> subdata;

};

PotData::PotData(bool progress_reports, int nr_subs, std::vector<double> *points, std::vector<double> *charges_coordinates, std::vector<double> *potential, pthread_mutex_t* mutex, int* progress_bar){
    m_nr_subs = nr_subs;
    m_nr_pnts = points->size();
    m_nr_ccos = charges_coordinates->size();
    m_nr_pots = potential->capacity(); //this has not yet been filled with values but has to be reserved already
    m_progress_bar = progress_bar;
    m_mut = mutex;
    m_progress_reports = progress_reports;
    //do a simple sanity check
    if (m_nr_pnts!=m_nr_pots*3){
        throw std::invalid_argument( "Points and potential do not have the same number of elements." );
    }
    //allocate data
    m_pnts = (double*)malloc(m_nr_pnts*sizeof(double));
    m_ccos = (double*)malloc(m_nr_ccos*sizeof(double));
    m_pots = (double*)malloc(m_nr_pots*sizeof(double));
    //copy data over
    //points:
    {
        double* temp = m_pnts;
        std::vector<double>::iterator it = points->begin();
        for (; it!=points->end(); ++it, ++temp){
            *temp = *it;
        } 
    }
    //charges_coordinates:
    {
        double* temp = m_ccos;
        std::vector<double>::iterator it = charges_coordinates->begin();
        for (; it!=charges_coordinates->end(); ++it, ++temp){
            *temp = *it;
        }
    }
    //create sub data
    subdata.reserve(m_nr_subs);
    const int min_nr_elements_per_sub = m_nr_pots/m_nr_subs;
    int too_many = m_nr_pots%m_nr_subs;
    int start = 0;
    for (int sub=0; sub<m_nr_subs; ++sub){
        const int here_elements = min_nr_elements_per_sub + ((too_many>0) ? 1 : 0);
        subdata.push_back(
                SubPotData(m_progress_reports, sub+1,
                           here_elements, m_nr_ccos, here_elements,
                           m_pnts+3*start,  m_ccos,    m_pots+start,
                           m_mut, m_progress_bar)
                );
        start += here_elements;
        too_many -=1;
    }
}

PotData::~PotData(){
    if (m_pnts!=NULL){
        free(m_pnts);
    }
    if (m_ccos!=NULL){
        free(m_ccos);
    }
    if (m_pots!=NULL){
        free(m_pots);
    }
}

PotData::PotData(){//parent constructor is called automatically if no arguments given
    m_nr_subs = 0;
}

SubPotData* PotData::GetSubPotData(int index){
    assert(index>=0 && index<m_nr_subs);
    return &(subdata[index]);
}

void PotData::TransferPotential(std::vector<double> *pot){
    double* temp = m_pots;
    for (int i=0; i<m_nr_pots; ++temp, ++i){
        pot->push_back(*temp);
    }
}
//END of class definitions

void* _potentialThread(void* data){
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED,NULL);
    struct timespec req = {0};
    req.tv_sec = 0;
    req.tv_nsec = 1L;
    SubPotData* dat = (SubPotData*) data;
    double *pnts = dat->GetPnts();
    double *ccos = dat->GetCCos();
    double *pots = dat->GetPots();
    const int nr_pnts = dat->GetNrPnts();
    const int nr_ccos = dat->GetNrCCos();
    const int nr_pots = dat->GetNrPots();
    const int sub_nr = dat->GetSubNr();
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
            nanosleep(&req, (struct timespec *)NULL);
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

void electrostatic_potential (bool progress_reports, int num_points, int num_charges, std::vector<double> points, std::vector<double> charges_coordinates, std::vector<double> *potential)
{
    int rc;
    int progress_bar = 0;
    if (progress_reports){
        printf("Starting: multi threaded computation of electrostatic potential.\n");
        fflush(stdout);
        rc = pthread_mutex_init(&mutex, NULL);
        if (rc){
            fprintf(stderr, "Failed to initialize mutex, will disable progress reports.\n");
            progress_reports = false;
        }
    }
    int num_threads = 1;
    {
        char const* tmp = getenv( "OMP_NUM_THREADS" );
        if ( tmp != NULL ) {
            num_threads = atoi(tmp);
        }
    }
    nr_threads = num_threads;
    threads = new pthread_t[num_threads];
    PotData *data;
    try
    {
        data = new PotData(progress_reports, num_threads, &points, &charges_coordinates, potential, &mutex, &progress_bar);
    }
    catch( const std::invalid_argument& e ) {
        throw;
    }
    signal(SIGINT, signal_callback_handler);
    for( int i=0; i < num_threads; ++i ){
        rc = pthread_create(threads+i, NULL, _potentialThread, (void *)data->GetSubPotData(i));
        assert(rc==0);
    }
    if ( progress_reports ){
        bool looping = true;
        int here_progress_bar = 0;
        fprintf(stdout, "Progress: %6.2f% | Total: %d/%d",0.0, 0, num_points);
        fflush(stdout);
        while (looping) {
            sleep(1);
            pthread_mutex_lock(&mutex);
            here_progress_bar = progress_bar;
            pthread_mutex_unlock(&mutex);
            printf("%c[2K\r", 27);
            if ( here_progress_bar >= num_points ){
                //print the correct information at the end
                fprintf(stdout,"Progress: %6.2f% | Total: %d/%d",100.0, num_points, num_points);
                looping = false;
            }
            else{
                //print intermittend information
                fprintf(stdout,"Progress: %6.2f% | Total: %d/%d",here_progress_bar*100.0/num_points, here_progress_bar, num_points);
            }
            fflush(stdout);
        }
        fprintf(stdout,"\n");
        fflush(stdout);
    }
    for( int i=0; i < num_threads; ++i ){
        pthread_join(*(threads+i), NULL);
    }
    data->TransferPotential(potential);
    delete[] threads;
    delete data;
    rc = pthread_mutex_trylock(&mutex);
    if (rc){
        fprintf(stderr, "The mutex could not be acquired, terminating abnormally.\n");
    }
    else{
        pthread_mutex_unlock(&mutex);
        pthread_mutex_destroy(&mutex);
    }
}
