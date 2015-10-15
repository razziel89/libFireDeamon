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
#ifndef H_PARALLEL_GENERIC_DEAMON_
#define H_PARALLEL_GENERIC_DEAMON_
#include <cstdlib>
#include <pthread.h>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <unistd.h>

//PG stands for ParallelGlobals
class PG{
    public:

        PG();
        ~PG();
        pthread_t* threads;
        pthread_mutex_t mutex;
        int nr_threads;
        int progress_bar;

};

extern PG* pg_global;

void signal_callback_handler(int signum);

//GPSubData stands for GenericParallelSubData
template <typename T>
class GPSubData{

    public:

        GPSubData();
        GPSubData(bool progress_reports, int sub_nr, std::vector<int> &len_data, std::vector<T*> &data, int len_output, T* output, pthread_mutex_t* mutex, int* progress_bar);
        ~GPSubData();
        T* GetData(int i);
        int GetNr(int i);
        T* GetDataOutput();
        int GetNrOutput();
        int GetSubNr();
        int* GetProgressBar();
        bool GetProgressReports();
        pthread_mutex_t* GetMutex();

    protected:

        std::vector<T*> m_data;
        std::vector<int> m_lengths;
        T* m_output;
        int m_len_output;
        int m_sub_nr;
        int* m_progress_bar;
        pthread_mutex_t* m_mut;
        bool m_progress_reports;

};

template <typename T>
GPSubData<T>::GPSubData(){
    m_data.reserve(0);
    m_lengths.reserve(0);
    m_output          = NULL;
    m_len_output      = 0;
    m_sub_nr          = 0;
    m_progress_bar    = NULL;
    m_mut             = NULL;
    m_progress_reports= false;
}

template <typename T>
GPSubData<T>::GPSubData(bool progress_reports, int sub_nr, std::vector<int> &lengths, std::vector<T*> &data, int len_output, T* output, pthread_mutex_t* mutex, int* progress_bar){
    m_progress_reports = progress_reports;
    m_mut              = mutex;
    m_progress_bar = progress_bar;
    m_sub_nr = sub_nr;
    m_output = output;
    m_len_output = len_output;
    m_lengths.reserve(lengths.size());
    for(std::vector<int>::const_iterator it = lengths.begin(); it != lengths.end(); ++it) {
        m_lengths.push_back(*it);
    }
    m_data.reserve(data.size());
    for(typename std::vector<T*>::const_iterator it = data.begin(); it != data.end(); ++it) {
        m_data.push_back(*it);
    }
}

template <typename T>
GPSubData<T>::~GPSubData(){}

template <typename T>
T* GPSubData<T>::GetData(int i){
    return m_data.at(i);
}

template <typename T>
int GPSubData<T>::GetNr(int i){
    return m_lengths.at(i);
}

template <typename T>
int GPSubData<T>::GetSubNr(){
    return m_sub_nr;
}

template <typename T>
int* GPSubData<T>::GetProgressBar(){
    return m_progress_bar;
}

template <typename T>
bool GPSubData<T>::GetProgressReports(){
    return m_progress_reports;
}

template <typename T>
pthread_mutex_t* GPSubData<T>::GetMutex(){
    return m_mut;
}

template <typename T>
T* GPSubData<T>::GetDataOutput(){
    return m_output;
}

template <typename T>
int GPSubData<T>::GetNrOutput(){
    return m_len_output;
}

template <typename T>
class GPData: public GPSubData<T> {

    public:

        GPData();
        GPData(bool progress_reports, int nr_subs, std::vector< std::vector<T> > &input, std::vector<T> *output, pthread_mutex_t* mutex, int* progress_bar, int split_index, int split_factor);
        ~GPData();
        GPSubData<T>* GetSubData(int i);
        void TransferOutput(bool empty_check = true);

    private:
        
        int m_nr_subs;
        int m_split_index;
        int m_split_factor;
        std::vector< GPSubData<T> > subdata;
        std::vector<T>* m_output_vector;
        using GPSubData<T>::m_data;
        using GPSubData<T>::m_lengths;
        using GPSubData<T>::m_output;
        using GPSubData<T>::m_len_output;
        using GPSubData<T>::m_sub_nr;
        using GPSubData<T>::m_progress_bar;
        using GPSubData<T>::m_mut;
        using GPSubData<T>::m_progress_reports;

};

template <typename T>
GPData<T>::GPData(bool progress_reports, int nr_subs, std::vector< std::vector<T> > &input, std::vector<T> *output, pthread_mutex_t* mutex, int* progress_bar, int split_index, int split_factor){
    m_nr_subs = nr_subs;
    m_lengths.reserve(input.size());
    m_data.reserve(input.size());
    for(typename std::vector< std::vector<T> >::const_iterator it = input.begin(); it != input.end(); ++it) {
        m_lengths.push_back(it->size());
    }
    m_len_output = output->capacity(); //this has not yet been filled with values but has to be reserved already
    m_progress_bar = progress_bar;
    m_mut = mutex;
    m_progress_reports = progress_reports;
    m_split_factor = split_factor;
    m_split_index = split_index;
    //do some sanity checking
    if (m_lengths[m_split_index]!=m_split_factor*m_len_output){
        throw std::invalid_argument( "The data to be split and the output data do not have the same number of elements considering the split factor." );
    }
    //allocate input data
    {
        std::vector<int>::iterator int_it = m_lengths.begin();
        for(; int_it != m_lengths.end(); ++int_it/*, ++T_it*/) {
            T* temp = (T*)malloc((*int_it)*sizeof(T));
            m_data.push_back(temp);
        }
    }
    //allocate output data
    m_output = (T*)malloc(m_len_output*sizeof(T));
    //remember the output vector
    m_output_vector = output;
    //copy input data over
    {
        typename std::vector<T*>::iterator T_it = m_data.begin();
        typename std::vector< std::vector<T> >::iterator input_it = input.begin();
        for(; T_it != m_data.end(); ++T_it, ++input_it) {
            T* temp = *T_it;
            typename std::vector<T>::iterator it = input_it->begin();
            for (; it!=input_it->end(); ++it, ++temp){
                *temp = *it;
            } 
        }
    }
    //create sub data
    subdata.reserve(m_nr_subs);
    const int min_nr_elements_per_sub = m_len_output/m_nr_subs;
    int too_many = m_len_output%m_nr_subs;
    int start = 0;
    for (int sub=0; sub<m_nr_subs; ++sub){
        const int here_elements = min_nr_elements_per_sub + ((too_many>0) ? 1 : 0);
        std::vector<int> here_len_data;
        here_len_data.reserve(m_lengths.size());
        typename std::vector<T*> here_data;
        here_data.reserve(m_data.size());
        {
            int count=0;
            std::vector<int>::iterator int_it = m_lengths.begin();
            typename std::vector<T*>::iterator T_it = m_data.begin();
            for(; int_it != m_lengths.end(); ++int_it, ++T_it, ++count) {
                if (count==m_split_index){
                    here_len_data.push_back(here_elements);
                    here_data.push_back((*T_it)+m_split_factor*start);
                }
                else{
                    here_len_data.push_back(*int_it);
                    here_data.push_back(*T_it);
                }
            }
        }
        subdata.push_back(
                GPSubData<T>(progress_reports, sub+1, //general information about thread
                             here_len_data, here_data, //input data
                             here_elements, m_output+start, //output data
                             m_mut, m_progress_bar) //information for parallelization
                );
        start += here_elements;
        too_many -=1;
    }
}

template <typename T>
GPData<T>::~GPData(){
    for(typename std::vector<T*>::iterator T_it = m_data.begin(); T_it != m_data.end(); ++T_it) {
        if (*T_it!=NULL){
            free(*T_it);
        }
    }
    if (m_output!=NULL){
        free(m_output);
    }
}

template <typename T>
GPData<T>::GPData(){//parent constructor is called automatically if no arguments given
    m_nr_subs = 0;
    m_split_index = 0;
    m_split_factor = 1;
    subdata.reserve(0);
}

template <typename T>
GPSubData<T>* GPData<T>::GetSubData(int i){
    assert(i>=0 && i<m_nr_subs);
    return &(subdata[i]);
}

template <typename T>
void GPData<T>::TransferOutput(bool empty_check){
    T* temp = m_output;
    if (empty_check && m_output_vector->size() > 0){
        throw std::length_error( "Output vector should be empty but it is not." );
    }
    for (int i=0; i<m_len_output; ++temp, ++i){
        m_output_vector->push_back(*temp);
    }
    //fprintf(stderr,"\nLengths: %d,   %d\n",m_output_vector->size(), m_len_output);fflush(stderr);
}
//END of class definitions

void init_parallel_generic(bool* progress_reports, PG* globals);

template <typename T>
void do_parallel_generic(void *(*thread_func)(void*), PG* globals, bool progress_reports, int nr_calcs, GPData<T>* data){
    if (pg_global){
        throw std::logic_error( "Global data for parallel computation already filled." );
    }
    else{
        pg_global = globals;
    }
    signal(SIGINT, signal_callback_handler);

    for( int i=0; i < globals->nr_threads; ++i ){
        int rc = pthread_create(globals->threads+i, NULL, thread_func, (void *)data->GetSubData(i));
        assert(rc==0);
    }
    if ( progress_reports ){
        bool looping = true;
        int here_progress_bar = 0;
        fprintf(stdout, "Progress: %6.2f%% | Total: %d/%d", 0.0 , 0, nr_calcs);
        fflush(stdout);
        while (looping) {
            sleep(1);
            pthread_mutex_lock(&(globals->mutex));
            here_progress_bar = globals->progress_bar;
            pthread_mutex_unlock(&(globals->mutex));
            fprintf(stdout,"%c[2K\r", 27);
            if ( here_progress_bar >= nr_calcs){
                fprintf(stdout,"Progress: %6.2f%% | Total: %d/%d", 100.0, nr_calcs, nr_calcs);
                looping = false;
            }
            else{
                fprintf(stdout,"Progress: %6.2f%% | Total: %d/%d",here_progress_bar*100.0/nr_calcs, here_progress_bar, nr_calcs);
            }
            fflush(stdout);
        }
        fprintf(stdout,"\n");
        fflush(stdout);
    }
    for( int i=0; i < globals->nr_threads; ++i ){
        pthread_join(*(globals->threads+i), NULL);
    }
}

void finalize_parallel_generic(bool progress_reports, PG* globals);

#endif //H_PARALLEL_GENERIC_DEAMON_
