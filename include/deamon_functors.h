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
#ifndef DEAMON_FUNCTORS_H
#define DEAMON_FUNCTORS_H

#include <vector>  //std::vector
#include <tuple>   //std::get
#include <cstring> //memcpy

//those functors can be used for the following things:

//add the sizes of each vector in a tuple to a vector
struct get_size_functor
{
    template<typename T, typename R>
    void operator () (T* t, R* r, int i)
    {
        r->push_back(t->size());
        (void) i; //suppress compiler warning
    }
};

//set each element in the tuple to NULL
struct set_to_NULL_functor
{
    template<typename T>
    void operator () (T* t, int i)
    {
        *t = NULL;
        (void) i; //suppress compiler warning
    }
};

//add the sizes in bytes of the data of each vector in a tuple to a vector
struct get_size_in_bytes_and_pointer_functor
{
    template<typename T, typename R>
    void operator () (std::vector<T>* t, R* r, int i)
    {
        T* data = t->data(); //WARNING: the error "void value not ignored as it ought to be" means you cannot use bool as a type, sorry
        r->push_back( std::make_tuple( t->size(), sizeof(T), (void*)data ) );
        (void) i; //suppress compiler warning
    }
};

//copy the data over
struct copy_functor_interlace
{
    int m_increment; //how many values belong together
    int m_nr_parts;  //in how many parts the data shall be split
    int m_nr_interlace; //which data stream to itnerlace (if at all)
    int m_interlace; //whether or not to interlace the data stream defined by m_nr_interlace

    copy_functor_interlace(int b, int s, int ni, bool i){
        m_increment    = b;
        m_nr_parts     = s;
        m_nr_interlace = ni;
        m_interlace    = i;
    }

    template<typename T, typename R>
    void operator () (T* t, R* r, int i)
    {

        //std::get<0>: number of elements in array
        //std::get<1>: size in bytes of data type
        //std::get<2>: pointers to the vectors data

        size_t total_size = std::get<0>((*r)[i]) * std::get<1>((*r)[i]);
        *t = (T)malloc(total_size);

        if (i == m_nr_interlace && m_interlace){
            size_t size = std::get<1>((*r)[i]) * m_increment;
            for (int part = 0; part < m_nr_parts; ++part){
                T start_pointer = (T)std::get<2>((*r)[i]);
                for (unsigned int dataset = 0; dataset < total_size/size; ++dataset){
                    if (dataset%part == 0){
                        memcpy( *t, start_pointer, size);
                        *t += m_increment;
                    }
                    start_pointer += m_increment;
                }
            }
        }
        else{
            memcpy( *t, std::get<2>((*r)[i]), total_size);
        }

    }

};

//free each element of the tuple
struct deallocate_functor
{
    template<typename T>
    void operator () (T* t, int i)
    {
        free(*t);
        //suppress compiler warnings about unused parameter
        //that I want to have turned off but these two parameters
        //are not needed here but needed in other functors
        (void) i;
    }
};

#endif //DEAMON_FUNCTORS_H
