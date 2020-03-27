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
/**
 * \file
 * \brief A header containing template classes and function definitions that allow to
 * perform parallelized computations.
 *
 * This is achieved by mapping a function to every data set in a std::vector (with
 * nigh-arbitrary template argument).  The type \a bool is not supported as either input
 * nor output type since std::vector<bool> is implemented as a bitfield and not as
 * simply avector of Boolean values.
 */
#ifndef H_PARALLEL_GENERIC_DEAMON_
#define H_PARALLEL_GENERIC_DEAMON_
#include <FireDeamon/core/deamon_functors.h>
#include <FireDeamon/core/iterate_over_tuple.h>
#include <assert.h>
#include <cstdlib>
#include <math.h>
#include <pthread.h>
#include <signal.h>
#include <stdexcept>
#include <stdio.h>
#include <time.h>
#include <tuple>
#include <unistd.h>
#include <vector>

#include <iostream>

// an alias template that is a tuple of vectors
template <typename... Ts> using tuple_of_vectors = std::tuple<std::vector<Ts>...>;

/**
 * \brief The class \a PG contains global information required for the parallelized
 * computation.
 *
 * The name stands for "ParallelGlobals". Parallelization is realized using multiple
 * threads via pthreads. A mutex (stands for "mutually exclusive") for manipulating
 * values in onjects of the class by all threads.  A progress bar is also provided to
 * allow printing progress reports.
 */
class PG {
public:
  //! \brief constructor
  PG();
  //! \brief destructor
  ~PG();
  /*! pointer to pthread_t - C-type array that allows for managing the threads
   * (contains thread handles)
   */
  pthread_t *threads;
  /*! pthread_mutex_t - a mutex that can be used to access data in a thread-safe way
   */
  pthread_mutex_t mutex;
  /*! int - the number of threads used for the parallel computation */
  int nr_threads;
  /*! int - a simple counter to estimate the progress of the compuation */
  int progress_bar;
};

extern PG *pg_global; // global instance

/**
 * \brief A function that is called whenever a signal is received (e.g., a keyboard
 * interrupt).
 *
 * Clean-up of data and thread-handles is also performed.
 */
void signal_callback_handler(int signum);

/**
 * \brief A templated class that contains all the data to be passed to single threads.
 *
 * Multiple instances of this class are aggregated in GPData. GPSubData stands for
 * "GenericParallelSubData". An arbitrary number of arguments can be passed to the
 * single threads.  Some bits for this class are taken from
 * http://stackoverflow.com/questions/27941661/generating-one-class-member-per-variadic-template-argument
 */
template <typename Tout, typename... Tins> class GPSubData {

public:
  //! \brief Default constructor
  GPSubData();
  /**
   * \brief Alternate constructor that allows to directly set most members.
   */
  GPSubData(
      /*! bool - whether or not a report on progress is desired */
      bool progress_reports,
      /*! int - a thread index (so that each threads knows its number in line) */
      int sub_nr,
      /*! std::vector<int> - a vector that contains the lengths of all elements in \a
       * data
       */
      std::vector<int> &len_data,
      /*! std::tuple<Tins*...> - a tuple aggregating all the data to be passed to the
       * threads. The data has to be in C-type array format.
       */
      std::tuple<Tins *...> &data,
      /*! int - the lengths of the output C-type array */
      int len_output,
      /*! pointer to Tout - this array will be filled with the output data */
      Tout *output,
      /*! pointer to pthread_mutex_t - the mutex to be used */
      pthread_mutex_t *mutex,
      /*! pointer to int - the counter used for reporting progress */
      int *progress_bar);
  //! \brief Default destructor
  ~GPSubData();
  /**
   * \brief A method to get the n-th set of input data.
   *
   * The number n is passed as a template argument.
   *
   * \return the n-th input C-type array
   */
  template <unsigned int N>
  typename std::tuple_element<N, std::tuple<Tins *...>>::type GetData();
  /**
   * \brief A method to get the number of entries in the n-th set of input data.
   *
   * The number n is passed as a template argument.
   *
   * \return the length of the n-th input C-type array
   */
  template <unsigned int N> int GetNr();
  //! \brief A method to get the C-type array for the output data
  //! \return C-type array for the output data
  Tout *GetDataOutput();
  //! \brief A method to get the length of the C-type array for the output data
  //! \return length of the C-type array for the output data
  int GetNrOutput();
  //! \brief Get the thread index.
  //! \return the thread index
  int GetSubNr();
  //! \brief Get the progress bar.
  //! \return a pointer to the int used to measure progress
  int *GetProgressBar();
  //! \brief Get \a progress_reports
  //! \return whether or not progress reports are desired
  bool GetProgressReports();
  //! \brief Get \a mutex
  //! \return a pointer to the mutex to be used
  pthread_mutex_t *GetMutex();

protected:
  /*! std::tuple<Tins*...> - the input data sets */
  std::tuple<Tins *...> m_data;
  /*! std::tuple<int> - the lenghts of the input data sets */
  std::vector<int> m_lengths;
  /*! pointer to Tout - C-type array for the output data */
  Tout *m_output;
  /*! int - length of the C-type array for the output dat */
  int m_len_output;
  /*! int - thread index */
  int m_sub_nr;
  /*! int - number of template arguments */
  int m_nr_types;
  /*! pointer to int - integer used to report progress */
  int *m_progress_bar;
  /*! pointer to pthread_mutex_t - mutex used for thread-safe access */
  pthread_mutex_t *m_mut;
  /*! bool - whether or not progress reports are desired */
  bool m_progress_reports;
};

template <typename Tout, typename... Tins> GPSubData<Tout, Tins...>::GPSubData() {
  m_lengths.reserve(0);
  m_nr_types = sizeof...(Tins);
  tuple_it::for_each_in_tuple(&m_data, set_to_NULL_functor());
  m_output = NULL;
  m_len_output = 0;
  m_sub_nr = 0;
  m_progress_bar = NULL;
  m_mut = NULL;
  m_progress_reports = false;
}

template <typename Tout, typename... Tins>
GPSubData<Tout, Tins...>::GPSubData(bool progress_reports, int sub_nr,
                                    std::vector<int> &lengths,
                                    std::tuple<Tins *...> &data, int len_output,
                                    Tout *output, pthread_mutex_t *mutex,
                                    int *progress_bar) {
  m_progress_reports = progress_reports;
  m_mut = mutex;
  m_progress_bar = progress_bar;
  m_sub_nr = sub_nr;
  m_output = output;
  m_len_output = len_output;
  m_lengths.reserve(lengths.size());
  for (std::vector<int>::const_iterator it = lengths.begin(); it != lengths.end();
       ++it) {
    m_lengths.push_back(*it);
  }
  m_data = data;
}

template <typename Tout, typename... Tins> GPSubData<Tout, Tins...>::~GPSubData() {}

template <typename Tout, typename... Tins>
template <unsigned int N>
typename std::tuple_element<N, std::tuple<Tins *...>>::type
GPSubData<Tout, Tins...>::GetData() {
  return std::get<N>(m_data);
}

template <typename Tout, typename... Tins>
template <unsigned int N>
int GPSubData<Tout, Tins...>::GetNr() {
  return m_lengths.at(N);
}

template <typename Tout, typename... Tins> int GPSubData<Tout, Tins...>::GetSubNr() {
  return m_sub_nr;
}

template <typename Tout, typename... Tins>
int *GPSubData<Tout, Tins...>::GetProgressBar() {
  return m_progress_bar;
}

template <typename Tout, typename... Tins>
bool GPSubData<Tout, Tins...>::GetProgressReports() {
  return m_progress_reports;
}

template <typename Tout, typename... Tins>
pthread_mutex_t *GPSubData<Tout, Tins...>::GetMutex() {
  return m_mut;
}

template <typename Tout, typename... Tins>
Tout *GPSubData<Tout, Tins...>::GetDataOutput() {
  return m_output;
}

template <typename Tout, typename... Tins> int GPSubData<Tout, Tins...>::GetNrOutput() {
  return m_len_output;
}

/**
 * \brief A templated class that contains all the data to be passed to all threads.
 *
 * This aggregates multiple instances of GPData and also spreads the data over all
 * threads. GPData stands for "GenericParallelData". An arbitrary number of arguments
 * can be passed to the single threads. Some bits for this class are taken from
 * http://stackoverflow.com/questions/27941661/generating-one-class-member-per-variadic-template-argument
 */
template <typename Tout, typename Tsplit, typename... Tins>
class GPData : public GPSubData<Tout, Tsplit, Tins...> {

public:
  //! \brief Default constructor
  GPData();
  /**
   * \brief Alternate constructor.
   *
   * Allows to set most members directly.
   */
  GPData(
      /*! bool - whether or not progress reports are desired */
      bool progress_reports,
      /*! int - the number of threads to use in parallel */
      int nr_subs,
      /**
       * \param input std::tuple<std::vector<Tsplit>,std::vector<Tins>...> - the input
       * data. The input data is given in the form of multiple objects of types
       * std::vector<Tins> and the vector whose content shall be spread over the
       * threads.
       */
      std::tuple<std::vector<Tsplit>, std::vector<Tins>...> &input,
      /*! pointer to std::vector<Tout> - the output data */
      std::vector<Tout> *output,
      /*! pointer to pthread_mutex_t - the mutex used to acces data thread-safely */
      pthread_mutex_t *mutex,
      /*! pointer to int - integer to be used to report progress */
      int *progress_bar,
      /*! int - the number of consecutive values in the vector (whose content is to be
       * spread over all threads) that shall remain together (e.g., would be 3 in the
       * case of Cartesian coordinates). Only used when interlacing.
       */
      int split_factor_in,
      /*! int - same as \a split_factor_in but for the output data */
      int split_factor_out,
      /*! bool - whether or not the input data shall be interlaced before being spread
       * over all threads. This might help to equalize loads.
       */
      bool interlace);
  //! \brief Default destructor
  ~GPData();
  //! \brief Get the i-th sub data (all data required for thread i in the form of an
  //! onject of type GPSubData)
  GPSubData<Tout, Tsplit, Tins...> *GetSubData(int i);
  //! \brief Transfer the output values from the C-type array to the std::vector<Tout>
  //! used for output
  void TransferOutput(bool empty_check = true);

private:
  /*! int - the number of threads to use in parallel */
  int m_nr_subs;
  /*! int - the number of consecutive values in the vector (whose content is to be
   * spread over all threads) that shall remain together (e.g., would be 3 in the case
   * of Cartesian coordinates). Only used when interlacing.
   */
  int m_split_factor_in;
  /*! int - same as \a m_split_factor_in but for the output data */
  int m_split_factor_out;
  /*! bool - whether or not the input data shall be interlaced before being spread
   * over all threads. This might help to equalize loads.
   */
  bool m_interlace;
  /*! std::vector<GPSubData<Tout,Tsplit,Tins...>> - data for all threads */
  std::vector<GPSubData<Tout, Tsplit, Tins...>> subdata;
  /*! pointer to std::vector<Tout> - will contain output data after calling \a
   * TransferOutput
   */
  std::vector<Tout> *m_output_vector;
  using GPSubData<Tout, Tsplit, Tins...>::m_data;
  using GPSubData<Tout, Tsplit, Tins...>::m_lengths;
  using GPSubData<Tout, Tsplit, Tins...>::m_output;
  using GPSubData<Tout, Tsplit, Tins...>::m_len_output;
  using GPSubData<Tout, Tsplit, Tins...>::m_sub_nr;
  using GPSubData<Tout, Tsplit, Tins...>::m_progress_bar;
  using GPSubData<Tout, Tsplit, Tins...>::m_mut;
  using GPSubData<Tout, Tsplit, Tins...>::m_progress_reports;
};

template <typename Tout, typename Tsplit, typename... Tins>
GPData<Tout, Tsplit, Tins...>::GPData(
    bool progress_reports, int nr_subs,
    std::tuple<std::vector<Tsplit>, std::vector<Tins>...> &input,
    std::vector<Tout> *output, pthread_mutex_t *mutex, int *progress_bar,
    int split_factor_in, int split_factor_out, bool interlace) {
  m_nr_subs = nr_subs > 0 ? nr_subs : 1;
  m_lengths.reserve(1 + sizeof...(Tins));
  // If there is only one thread, do not waste time with interlacing and
  // deinterlacing since that is not necessary
  m_interlace = m_nr_subs > 1 ? interlace : false;
  tuple_it::for_each_in_tuple_vector(&input, &m_lengths, get_size_functor());
  // The following has not yet been filled with values but has to be reserved already
  m_len_output = output->capacity();
  m_progress_bar = progress_bar;
  m_mut = mutex;
  m_progress_reports = progress_reports;
  m_split_factor_in = split_factor_in;
  m_split_factor_out = split_factor_out;
  // do some sanity checking
  if (m_len_output % m_split_factor_out != 0) {
    throw std::invalid_argument("The number of elements in the output data is not "
                                "divisible by the output split factor.");
  }
  if (m_lengths[0] % m_split_factor_in != 0) {
    throw std::invalid_argument("The number of elements in the input data (split "
                                "coloumn) is not divisible by the input split factor.");
  }
  if (m_lengths[0] / m_split_factor_in != m_len_output / m_split_factor_out) {
    throw std::invalid_argument(
        "The data to be split and the output data do not have the same number of "
        "elements considering the split factor.");
  }
  // allocate output data
  m_output = (Tout *)malloc(m_len_output * sizeof(Tout));
  // remember the output vector
  m_output_vector = output;
  // copy input data over and allocate beforehand
  {
    std::vector<std::tuple<unsigned int, size_t, void *>> sizes_pointers_vec;
    sizes_pointers_vec.reserve(1 + sizeof...(Tins));

    tuple_it::for_each_in_tuple_vector(
        &input, &sizes_pointers_vec, get_size_in_bytes_and_pointer_functor());

    tuple_it::for_each_in_tuple_vector(
        &m_data,
        &sizes_pointers_vec,
        copy_functor_interlace(m_split_factor_in, m_nr_subs, 0, m_interlace));
  }
  // create sub data
  subdata.reserve(m_nr_subs);
  const int min_nr_elements_per_sub = (m_len_output / m_split_factor_out) / m_nr_subs;
  int too_many = (m_len_output / m_split_factor_out) % m_nr_subs;
  int start = 0;
  for (int sub = 0; sub < m_nr_subs; ++sub) {
    const int here_elements = min_nr_elements_per_sub + ((too_many > 0) ? 1 : 0);
    std::vector<int> here_len_data;
    here_len_data.reserve(m_lengths.size());
    // copy assignment
    typename std::tuple<Tsplit *, Tins *...> here_data = m_data;
    {
      int count = 0;
      std::vector<int>::iterator int_it = m_lengths.begin();
      // typename std::vector<Tin*>::iterator T_it = m_data.begin();
      for (; int_it != m_lengths.end(); ++int_it /*, ++T_it*/, ++count) {
        if (count == 0) {
          here_len_data.push_back(here_elements);
          // here_data.push_back((*T_it)+m_split_factor_in*start);
          // push the vector for the input stream which is being split over
          // threads forward by the necessary number of steps
          std::get<0>(here_data) += m_split_factor_in * start;
        } else {
          here_len_data.push_back(*int_it);
          // here_data.push_back(*T_it);
        }
      }
    }
    subdata.push_back(GPSubData<Tout, Tsplit, Tins...>(
        // general information about thread
        progress_reports,
        sub + 1,
        // input data
        here_len_data,
        here_data,
        // output data
        here_elements,
        m_output + m_split_factor_out * start,
        // information for parallelization
        m_mut,
        m_progress_bar));
    start += here_elements;
    too_many -= 1;
  }
}

template <typename Tout, typename Tsplit, typename... Tins>
GPData<Tout, Tsplit, Tins...>::~GPData() {
  tuple_it::for_each_in_tuple(&m_data, deallocate_functor());
  if (m_output != NULL) {
    free(m_output);
  }
}

template <typename Tout, typename Tsplit, typename... Tins>
GPData<Tout, Tsplit, Tins...>::GPData() {
  // parent constructor is called automatically if no arguments given
  m_nr_subs = 0;
  m_split_factor_in = 1;
  m_split_factor_out = 1;
  m_interlace = false;
  subdata.reserve(0);
  m_output_vector.reserve(0);
}

template <typename Tout, typename Tsplit, typename... Tins>
GPSubData<Tout, Tsplit, Tins...> *GPData<Tout, Tsplit, Tins...>::GetSubData(int i) {
  assert(i >= 0 && i < m_nr_subs);
  return &(subdata[i]);
}

template <typename Tout, typename Tsplit, typename... Tins>
void GPData<Tout, Tsplit, Tins...>::TransferOutput(bool empty_check) {
  if (empty_check && m_output_vector->size() > 0) {
    throw std::length_error("Output vector should be empty but it is not.");
  }
  if (m_interlace) { // deinterlace the output data
    Tout **temp = (Tout **)malloc(m_nr_subs * sizeof(Tout *));
    const int min_nr_elements_per_sub = (m_len_output / m_split_factor_out) / m_nr_subs;
    int too_many = (m_len_output / m_split_factor_out) % m_nr_subs;
    int start = 0;
    for (int sub = 0; sub < m_nr_subs; ++sub) {
      temp[sub] = m_output + (m_split_factor_out * start);
      const int here_elements = min_nr_elements_per_sub + ((too_many > 0) ? 1 : 0);
      start += here_elements;
      too_many -= 1;
    }
    for (int i = 0; i < m_len_output / m_split_factor_out; ++i) {
      int modulus = i % m_nr_subs;
      for (int j = 0; j < m_split_factor_out; ++j) {
        m_output_vector->push_back(*(temp[modulus]));
        ++(temp[modulus]);
      }
    }
    free(temp);
  } else {
    Tout *temp = m_output;
    for (int i = 0; i < m_len_output; ++temp, ++i) {
      m_output_vector->push_back(*temp);
    }
  }
}
// END of class definitions

/**
 * \brief initialize the global data structure (that is used for signal handling and
 * reporting progress)
 */
void init_parallel_generic(bool *progress_reports, PG *globals);

/**
 * \brief Perform a parallelized compuation.
 */
template <typename... Ts>
void do_parallel_generic(
    /*! void *(*thread_func)(void*) - function pointer. This function is mapped to the
     * data.
     */
    void *(*thread_func)(void *),
    /*! pointer to PG - global data that is, e.g., used for treating keyboard interrupts
     */
    PG *globals,
    /*! bool - whether or not progress reports are desired */
    bool progress_reports,
    /*! int - how many computations shall be performed ,i.e., maximum counter for
     * progress reports
     */
    int nr_calcs,
    /*! pointer to GPData<Ts...> - the data structure that contains all the data */
    GPData<Ts...> *data) {
  if (pg_global) {
    throw std::logic_error("Global data for parallel computation already filled.");
  } else {
    pg_global = globals;
  }
  signal(SIGINT, signal_callback_handler);

  for (int i = 0; i < globals->nr_threads; ++i) {
    int rc = pthread_create(
        globals->threads + i, NULL, thread_func, (void *)data->GetSubData(i));
    assert(rc == 0);
  }
  char prog_char = '|';
  if (progress_reports) {
    bool looping = true;
    int here_progress_bar = 0;
    fprintf(stdout, "Progress: %6.2f%% | Total: %d/%d %c", 0.0, 0, nr_calcs, prog_char);
    fflush(stdout);
    while (looping) {
      sleep(1);
      switch (prog_char) {
      case '|':
        prog_char = '/';
        break;
      case '/':
        prog_char = '-';
        break;
      case '-':
        prog_char = '\\';
        break;
      case '\\':
        prog_char = '|';
        break;
      }
      pthread_mutex_lock(&(globals->mutex));
      here_progress_bar = globals->progress_bar;
      pthread_mutex_unlock(&(globals->mutex));
      fprintf(stdout, "%c[2K\r", 27);
      if (here_progress_bar >= nr_calcs) {
        fprintf(
            stdout, "Progress: %6.2f%% | Total: %d/%d  ", 100.0, nr_calcs, nr_calcs);
        looping = false;
      } else {
        fprintf(stdout,
                "Progress: %6.2f%% | Total: %d/%d %c",
                here_progress_bar * 100.0 / nr_calcs,
                here_progress_bar,
                nr_calcs,
                prog_char);
      }
      fflush(stdout);
    }
    fprintf(stdout, "\n");
    fflush(stdout);
  }
  for (int i = 0; i < globals->nr_threads; ++i) {
    pthread_join(*(globals->threads + i), NULL);
  }
}

/**
 * \brief finalize everything after the parallel computation. This also transfers output
 * data properly.
 */
void finalize_parallel_generic(bool progress_reports, PG *globals);

#endif // H_PARALLEL_GENERIC_DEAMON_
