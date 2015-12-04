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
#include <tuple>
#include <stdexcept>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <unistd.h>
#include <iterate_over_tuple.h>
#include <deamon_functors.h>

//an alias template that is a tuple of vectors
template <typename... Ts>
using tuple_of_vectors = std::tuple<std::vector<Ts>...>;

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
//some bits taken from http://stackoverflow.com/questions/27941661/generating-one-class-member-per-variadic-template-argument
template <typename Tout, typename... Tins>
class GPSubData{

    public:

        GPSubData();
        GPSubData(bool progress_reports, int sub_nr, std::vector<int> &len_data, std::tuple<Tins*...> &data, int len_output, Tout* output, pthread_mutex_t* mutex, int* progress_bar);
        ~GPSubData();
        template <unsigned int N>
        typename std::tuple_element<N,std::tuple<Tins*...>>::type GetData();
        template <unsigned int N>
        int GetNr();
        Tout* GetDataOutput();
        int GetNrOutput();
        int GetSubNr();
        int* GetProgressBar();
        bool GetProgressReports();
        pthread_mutex_t* GetMutex();

    protected:

        std::tuple<Tins*...> m_data;
        std::vector<int> m_lengths;
        Tout* m_output;
        int m_len_output;
        int m_sub_nr;
        int m_nr_types;
        int* m_progress_bar;
        pthread_mutex_t* m_mut;
        bool m_progress_reports;

};

template <typename Tout, typename... Tins>
GPSubData<Tout,Tins...>::GPSubData(){
    m_lengths.reserve(0);
    m_nr_types        = sizeof...(Tins);
    for_each_in_tuple(&m_data,set_to_NULL_functor());
    m_output          = NULL;
    m_len_output      = 0;
    m_sub_nr          = 0;
    m_progress_bar    = NULL;
    m_mut             = NULL;
    m_progress_reports= false;
}

template <typename Tout, typename... Tins>
GPSubData<Tout,Tins...>::GPSubData(bool progress_reports, int sub_nr, std::vector<int> &lengths, std::tuple<Tins*...> &data, int len_output, Tout* output, pthread_mutex_t* mutex, int* progress_bar){
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
    m_data = data;
}

template <typename Tout, typename... Tins>
GPSubData<Tout,Tins...>::~GPSubData(){}

template <typename Tout, typename... Tins>
    template <unsigned int N>
typename std::tuple_element<N,std::tuple<Tins*...>>::type GPSubData<Tout,Tins...>::GetData(){
    return std::get<N>(m_data);
}

template <typename Tout, typename... Tins>
    template <unsigned int N>
int GPSubData<Tout,Tins...>::GetNr(){
    return m_lengths.at(N);
}

template <typename Tout, typename... Tins>
int GPSubData<Tout,Tins...>::GetSubNr(){
    return m_sub_nr;
}

template <typename Tout, typename... Tins>
int* GPSubData<Tout,Tins...>::GetProgressBar(){
    return m_progress_bar;
}

template <typename Tout, typename... Tins>
bool GPSubData<Tout,Tins...>::GetProgressReports(){
    return m_progress_reports;
}

template <typename Tout, typename... Tins>
pthread_mutex_t* GPSubData<Tout,Tins...>::GetMutex(){
    return m_mut;
}

template <typename Tout, typename... Tins>
Tout* GPSubData<Tout,Tins...>::GetDataOutput(){
    return m_output;
}

template <typename Tout, typename... Tins>
int GPSubData<Tout,Tins...>::GetNrOutput(){
    return m_len_output;
}

template <typename Tout, typename Tsplit, typename... Tins>
class GPData: public GPSubData<Tout, Tsplit, Tins...> {

    public:

        GPData();
        GPData(bool progress_reports, int nr_subs, std::tuple< std::vector<Tsplit>,std::vector<Tins>...> &input, std::vector<Tout> *output, pthread_mutex_t* mutex, int* progress_bar, int split_factor_in, int split_factor_out, bool interlace);
        ~GPData();
        GPSubData<Tout,Tsplit,Tins...>* GetSubData(int i);
        void TransferOutput(bool empty_check = true);

    private:
        
        int m_nr_subs;
        int m_split_factor_in;
        int m_split_factor_out;
        bool m_interlace;
        std::vector< GPSubData<Tout,Tsplit,Tins...> > subdata;
        std::vector<Tout>* m_output_vector;
        using GPSubData<Tout,Tsplit,Tins...>::m_data;
        using GPSubData<Tout,Tsplit,Tins...>::m_lengths;
        using GPSubData<Tout,Tsplit,Tins...>::m_output;
        using GPSubData<Tout,Tsplit,Tins...>::m_len_output;
        using GPSubData<Tout,Tsplit,Tins...>::m_sub_nr;
        using GPSubData<Tout,Tsplit,Tins...>::m_progress_bar;
        using GPSubData<Tout,Tsplit,Tins...>::m_mut;
        using GPSubData<Tout,Tsplit,Tins...>::m_progress_reports;

};

template <typename Tout, typename Tsplit, typename... Tins>
GPData<Tout,Tsplit,Tins...>::GPData(bool progress_reports, int nr_subs, std::tuple< std::vector<Tsplit>, std::vector<Tins>...> &input, std::vector<Tout> *output, pthread_mutex_t* mutex, int* progress_bar, int split_factor_in, int split_factor_out, bool interlace){
    m_nr_subs = nr_subs>0 ? nr_subs : 1;
    m_lengths.reserve(1 + sizeof...(Tins));
    //if there is only one thread, do not waste time with interlacing and
    //deinterlacing since that is not necessary
    m_interlace = m_nr_subs>1 ? interlace : false;
    for_each_in_tuple_vector(&input, &m_lengths, get_size_functor());
    m_len_output = output->capacity(); //this has not yet been filled with values but has to be reserved already
    m_progress_bar = progress_bar;
    m_mut = mutex;
    m_progress_reports = progress_reports;
    m_split_factor_in  = split_factor_in;
    m_split_factor_out = split_factor_out;
    //do some sanity checking
    if (m_len_output%m_split_factor_out != 0){
        throw std::invalid_argument( "The number of elements in the output data is not divisible by the output split factor." );
    }
    if (m_lengths[0]%m_split_factor_in != 0){
        throw std::invalid_argument( "The number of elements in the input data (split coloumn) is not divisible by the input split factor." );
    }
    if (m_lengths[0]/m_split_factor_in != m_len_output/m_split_factor_out){
        throw std::invalid_argument( "The data to be split and the output data do not have the same number of elements considering the split factor." );
    }
    //allocate output data
    m_output = (Tout*)malloc(m_len_output*sizeof(Tout));
    //remember the output vector
    m_output_vector = output;
    //copy input data over and allocate beforehand
    {
        std::vector<std::tuple<unsigned int,size_t,void*>> sizes_pointers_vec;
        sizes_pointers_vec.reserve(1 + sizeof...(Tins));

        for_each_in_tuple_vector(&input, &sizes_pointers_vec, get_size_in_bytes_and_pointer_functor());

        for_each_in_tuple_vector(&m_data, &sizes_pointers_vec, copy_functor_interlace(m_split_factor_in,m_nr_subs,m_interlace,0));
    }
    //create sub data
    subdata.reserve(m_nr_subs);
    const int min_nr_elements_per_sub = (m_len_output/m_split_factor_out)/m_nr_subs;
    int too_many = (m_len_output/m_split_factor_out)%m_nr_subs;
    int start = 0;
    for (int sub=0; sub<m_nr_subs; ++sub){
        const int here_elements = min_nr_elements_per_sub + ((too_many>0) ? 1 : 0);
        std::vector<int> here_len_data;
        here_len_data.reserve(m_lengths.size());
        //copy assignment
        typename std::tuple<Tsplit*,Tins*...> here_data = m_data;
        {
            int count=0;
            std::vector<int>::iterator int_it = m_lengths.begin();
            //typename std::vector<Tin*>::iterator T_it = m_data.begin();
            for(; int_it != m_lengths.end(); ++int_it/*, ++T_it*/, ++count) {
                if (count==0){
                    here_len_data.push_back(here_elements);
                    //here_data.push_back((*T_it)+m_split_factor_in*start);
                    //push the vector for the input stream which is being split over
                    //threads forward by the necessary number of steps
                    std::get<0>(here_data) += m_split_factor_in*start;
                }
                else{
                    here_len_data.push_back(*int_it);
                    //here_data.push_back(*T_it);
                }
            }
        }
        subdata.push_back(
                GPSubData<Tout,Tsplit,Tins...>(
                                    progress_reports, sub+1, //general information about thread
                                    here_len_data, here_data, //input data
                                    here_elements, m_output+m_split_factor_out*start, //output data
                                    m_mut, m_progress_bar) //information for parallelization
                );
        start += here_elements;
        too_many -=1;
    }
}

template <typename Tout, typename Tsplit, typename... Tins>
GPData<Tout,Tsplit,Tins...>::~GPData(){
    for_each_in_tuple(&m_data,deallocate_functor());
    if (m_output!=NULL){
        free(m_output);
    }
}

template <typename Tout, typename Tsplit, typename... Tins>
GPData<Tout,Tsplit,Tins...>::GPData(){//parent constructor is called automatically if no arguments given
    m_nr_subs = 0;
    m_split_factor_in = 1;
    m_split_factor_out = 1;
    m_interlace = false;
    subdata.reserve(0);
    m_output_vector.reserve(0);
}

template <typename Tout, typename Tsplit, typename... Tins>
GPSubData<Tout,Tsplit,Tins...>* GPData<Tout,Tsplit,Tins...>::GetSubData(int i){
    assert(i>=0 && i<m_nr_subs);
    return &(subdata[i]);
}

template <typename Tout, typename Tsplit, typename... Tins>
void GPData<Tout,Tsplit,Tins...>::TransferOutput(bool empty_check){
    if (empty_check && m_output_vector->size() > 0){
        throw std::length_error( "Output vector should be empty but it is not." );
    }
    if (m_interlace){//deinterlace the output data
        Tout** temp = (Tout**)malloc(m_nr_subs*sizeof(Tout*));
        const int min_nr_elements_per_sub = (m_len_output/m_split_factor_out)/m_nr_subs;
        int too_many = (m_len_output/m_split_factor_out)%m_nr_subs;
        int start = 0;
        for (int sub=0; sub<m_nr_subs; ++sub){
            temp[sub] = m_output+(m_split_factor_out*start);
            const int here_elements = min_nr_elements_per_sub + ((too_many>0) ? 1 : 0);
            start += here_elements;
            too_many -=1;
        }
        for (int i=0; i<m_len_output/m_split_factor_out; ++i){
            int modulus = i%m_nr_subs;
            for (int j=0; j<m_split_factor_out; ++j){
                m_output_vector->push_back(*(temp[modulus]));
                ++(temp[modulus]);
            }
        }
        free(temp);
    }
    else{
        Tout* temp = m_output;
        for (int i=0; i<m_len_output; ++temp, ++i){
            m_output_vector->push_back(*temp);
        }
    }
}
//END of class definitions

void init_parallel_generic(bool* progress_reports, PG* globals);

template <typename... Ts>
void do_parallel_generic(void *(*thread_func)(void*), PG* globals, bool progress_reports, int nr_calcs, GPData<Ts...>* data){
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
