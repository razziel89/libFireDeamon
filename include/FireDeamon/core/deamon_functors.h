/***********
This file is part of libFireDeamon.

Copyright (C) 2016 by Torsten Sachse

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

#include <cstring> //memcpy
#include <tuple>   //std::get
#include <vector>  //std::vector
/**
 * \file
 * \brief A header that containes some functors that allow to do some things
 * for each entry in a tuple.
 *
 * They are used in conjunction with iterate_over_tuple.h to
 * do that.
 */
/// \name Functors
/// \{
struct get_size_functor;
struct set_to_NULL_functor;
struct get_size_in_bytes_and_pointer_functor;
struct copy_functor_interlace;
struct deallocate_functor;
/// \}

//! \brief Add the sizes of a vector to a vector
struct get_size_functor {
  /**
   * \brief Operator that performs the operation
   * \param t std::vector<T>* - pointer to the vector whose length shall be added to a
   * vector
   * \param i int - helper parameter that allows for looping over each element in a
   * tuple
   */
  template <typename T> void operator()(std::vector<T> *t, std::vector<int> *r, int i) {
    r->push_back(t->size());
    (void)i; // suppress compiler warning
  }
};

//! \brief Set a pointer to NULL
struct set_to_NULL_functor {
  /**
   * \brief Operator that performs the operation
   * \param t T** - pointer to a pointer that shall be set to NULL
   * \param i int - helper parameter that allows for looping over each element in a
   * tuple
   */
  template <typename T> void operator()(T **t, int i) {
    *t = NULL;
    (void)i; // suppress compiler warning
  }
};

//! \brief Create a tuple that contains information about a vector and append that tuple
//!  to a vector.
//!
//! The information contained in the tuple that is creates is as follows:
//! -# number of elements in the vector
//! -# size in bytes of data type
//! -# pointers to the vector's data
struct get_size_in_bytes_and_pointer_functor {
  /**
   * \brief Operator that performs the operation
   * \param t std::vector<T>* - pointer to a vecotor whose information shall be
   * extracted
   * \param r std::vector<std::tuple<unsigned int,size_t,void*>>* - pointer to the
   * vector to which to append the tuple
   * \param i int - helper parameter that allows for looping over each element in a
   * tuple
   */
  template <typename T>
  void operator()(std::vector<T> *t,
                  std::vector<std::tuple<unsigned int, size_t, void *>> *r, int i) {
    unsigned int temp_vect_size = t->size();
    T *data = t->data(); // WARNING: the error "void value not ignored as it ought
                         // to be" means you cannot use bool as a type, sorry
    size_t T_size = sizeof(T);
    void *temp_void_ptr = (void *)data;
    std::tuple<unsigned int, size_t, void *> temp_tuple =
        std::make_tuple(temp_vect_size, T_size, temp_void_ptr);
    r->push_back(temp_tuple);
    (void)i; // suppress compiler warning
  }
};

//! \brief Copy the data in a vector over to a number of C-type arrays each (supports
//!  interlacing)
//!
//! Data can be grouped together, meaing that it is possible to keep a set of data
//! together even when interlacing the data during the copy process. Interlacing the
//! data can help to balance the load when performing computations.
struct copy_functor_interlace {
  /*! Size of the group that belongs together (important when interlacing) */
  unsigned int m_increment;
  /*! \brief In how many parts the data shall be split, i.e., how many threads will be
   * used for parallel computations
   */
  unsigned int m_nr_parts;
  /*! Index of data stream to interlace (if at all) */
  int m_nr_interlace;
  /*! Whether or not to interlace the data stream defined by m_nr_interlace */
  int m_interlace;

  /**
   * \brief Constructor for the functor
   * \param b int - Size of the group that belongs together (important when
   * interlacing)
   * \param s int - In how many parts the data shall be split, i.e., how many threads
   * will perform a computation simultaneously
   * \param ni int - Index of data stream to interlace (if at all)
   * \param i bool - Whether or not to interlace the data stream defined by
   * m_nr_interlace
   */
  copy_functor_interlace(int b, int s, int ni, bool i) {
    m_increment = (unsigned int)b;
    m_nr_parts = (unsigned int)s;
    m_nr_interlace = ni;
    m_interlace = i;
  }

  /**
   * \brief Operator that performs the operation.
   * \param t T** - pointer to C-type array to which the data shall be copied
   * \param r std::vector<std::tuple<unsigned int,size_t,void*>>* - a vector
   * containing the information about the data that is to be copied over (generated by
   * \a get_size_in_bytes_and_pointer_functor)
   * \param i int - helper parameter that allows for looping over each element in a
   * tuple (this is also the index for the data taken from \a r)
   */
  template <typename T>
  void operator()(T **t, std::vector<std::tuple<unsigned int, size_t, void *>> *r,
                  int i) {

    // std::get<0>: number of elements in array
    // std::get<1>: size in bytes of data type
    // std::get<2>: pointers to the vectors data

    unsigned int nr_elements = std::get<0>((*r)[i]);
    size_t element_size = std::get<1>((*r)[i]);
    size_t total_size = nr_elements * element_size;
    *t = (T *)malloc(total_size);

    if (i == m_nr_interlace && m_interlace) {
      size_t size_dataset = element_size * m_increment;
      T *dest_pointer = *t;
      for (unsigned int part = 0; part < m_nr_parts; ++part) {
        T *src_pointer = (T *)std::get<2>((*r)[i]);
        for (unsigned int src_dataset = 0; src_dataset < nr_elements / m_increment;
             ++src_dataset, src_pointer += m_increment) {
          if (src_dataset % m_nr_parts == part) {
            memcpy(dest_pointer, src_pointer, size_dataset);
            dest_pointer += m_increment;
          }
        }
      }
      // the commented out code block below performs interlacing without
      ////unsigned int dest_dataset = 0;
      // T dest_pointer = *t;
      // for (unsigned int part = 0; part < m_nr_parts; ++part){
      //    T src_pointer  = (T)std::get<2>((*r)[i]);
      //    for (unsigned int src_dataset = 0; src_dataset < nr_elements; ++src_dataset,
      //    ++src_pointer){
      //        if ((src_dataset/m_increment)%m_nr_parts == part){
      //            *dest_pointer = *src_pointer;
      //            ++dest_pointer;
      //        }
      //    }
      //}
    } else {
      memcpy(*t, std::get<2>((*r)[i]), total_size);
    }
  }
};

//! \brief Free each element in the tuple
struct deallocate_functor {
  /**
   * \brief Operator that performs the operation
   * \param t T** - pointer to a pointer that shall be freed
   * \param i int - helper parameter that allows for looping over each element in a
   * tuple
   */
  template <typename T> void operator()(T **t, int i) {
    free(*t);
    // suppress compiler warnings about unused parameter that I want to have turned
    // off but these two parameters are not needed here but needed in other functors
    (void)i; // suppress compiler warning
  }
};
#endif // DEAMON_FUNCTORS_H
