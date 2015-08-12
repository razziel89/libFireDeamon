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

//BEGINNING of class definitions
//this class holds all the data that will be worked on
//I have chosen to use a class in order to be able to pass a single argument to the thread functions
class SubPotData {

    public:

        SubPotData();
        SubPotData(int progress, int sub_nr, int l_pnts, int l_ccos, int l_pots, double* points, double* charges_coordinates, double* potential);
        ~SubPotData();
        double* GetPnts();
        double* GetCCos();
        double* GetPots();
        int GetNrPnts();
        int GetNrCCos();
        int GetNrPots();
        int GetSubNr();
        int GetProgress();

    protected:

        double *m_pnts;
        double *m_ccos;
        double *m_pots;
        int m_nr_pnts;
        int m_nr_ccos;
        int m_nr_pots;
        int m_sub_nr;
        int m_progress;

};

SubPotData::SubPotData(){
    m_pnts     = NULL;
    m_ccos     = NULL;
    m_pots     = NULL;
    m_nr_pnts  = 0;
    m_nr_ccos  = 0;
    m_nr_pots  = 0;
    m_sub_nr   = 0;
    m_progress = 0;
}

SubPotData::SubPotData(int progress, int sub_nr, int l_pnts, int l_ccos, int l_pots, double* points, double* charges_coordinates, double* potential){
    m_nr_pnts = l_pnts;
    m_nr_ccos = l_ccos;
    m_nr_pots = l_pots;
    m_sub_nr = sub_nr;
    m_progress = progress;
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
int SubPotData::GetProgress(){
    return m_progress;
}

class PotData : public SubPotData {

    public:

        PotData();
        PotData(int progress, int nr_subs, std::vector<double> *points, std::vector<double> *charges_coordinates, std::vector<double> *potential);
        ~PotData();
        SubPotData* GetSubPotData(int index);
        void TransferPotential(std::vector<double> *pot);

    private:
        
        int m_nr_subs;
        std::vector<SubPotData> subdata;

};

PotData::PotData(int progress, int nr_subs, std::vector<double> *points, std::vector<double> *charges_coordinates, std::vector<double> *potential){
    m_nr_subs = nr_subs;
    m_nr_pnts = points->size();
    m_nr_ccos = charges_coordinates->size();
    m_nr_pots = potential->capacity(); //this has not yet been filled with values but has to be reserved already
    m_progress = progress;
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
                SubPotData(m_progress, sub+1,
                           here_elements, m_nr_ccos, here_elements,
                           m_pnts+3*start,  m_ccos,    m_pots+start)
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
    SubPotData* dat = (SubPotData*) data;
    double *pnts = dat->GetPnts();
    double *ccos = dat->GetCCos();
    double *pots = dat->GetPots();
    const int nr_pnts = dat->GetNrPnts();
    const int nr_ccos = dat->GetNrCCos();
    const int nr_pots = dat->GetNrPots();
    const int progress = dat->GetProgress();
    const int sub_nr = dat->GetSubNr();

    double* p = pots;
    double* pc = pnts;
    if ( progress > 0 ){
        printf("Thread: %2d | Progress: %5.2f% | Nr. computations per progress report: %li\n",sub_nr,0.0,nr_ccos*progress);
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
            printf("Thread: %2d | Progress: %5.2f%\n",sub_nr,100.0*ip/nr_pots);
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
        }
    }

    pthread_exit(NULL);
}

void electrostatic_potential (bool progress_reports, int num_points, int num_charges, std::vector<double> points, std::vector<double> charges_coordinates, std::vector<double> *potential)
{
    int progress = 0;
    if (progress_reports){
        progress = num_points/250;
        printf("Starting: multi threaded computation of electrostatic potential.\n");
        fflush(stdout);
    }
    int num_threads = 1;
    {
        char const* tmp = getenv( "OMP_NUM_THREADS" );
        if ( tmp != NULL ) {
            num_threads = atoi(tmp);
        }
    }
    PotData *data;
    try
    {
        data = new PotData(progress, num_threads, &points, &charges_coordinates, potential);
    }
    catch( const std::invalid_argument& e ) {
        throw;
    }
    pthread_t threads[num_threads];
    int rc;
    for( int i=0; i < num_threads; ++i ){
        rc = pthread_create(&threads[i], NULL, _potentialThread, (void *)data->GetSubPotData(i));
        assert(rc==0);
    }
    for( int i=0; i < num_threads; ++i ){
        pthread_join(threads[i], NULL);
    }
    data->TransferPotential(potential);
    delete data;
}
