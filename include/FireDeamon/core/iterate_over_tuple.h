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
 * \brief Header file aiding in executing code for every entry in a tuple.
 */
#ifndef ITERATE_OVER_TUPLE_DEAMON_H
#define ITERATE_OVER_TUPLE_DEAMON_H

#include <tuple>

// taken from
// http://stackoverflow.com/questions/16387354/template-tuple-calling-a-function-on-each-element
/**
 * \brief namespace containing templates that can be used to perform actions for every
 * entry in a tuple.
 *
 * The "iteration" is no actual iteration as tuples are objects whose lengths have to be
 * fully known at compile time. Access functions also have to be known at compile time.
 * Hence, templating is used to create a sequence 1..N where N is the length of the
 * tuple over which to "iterate".
 */
namespace tuple_it {
//! \brief generate a sequence of numbers
template <int... Is> struct seq {};

//! \brief recursively generate a sequence of numbers and keep them in the template
//! information
template <int N, int... Is> struct gen_seq : gen_seq<N - 1, N - 1, Is...> {};

//! \brief the struct that is the end of the recursion
template <int... Is> struct gen_seq<0, Is...> : seq<Is...> {};

/**
 * \brief Evaluate the functor for each element of the tuple. Not to be called
 * directly.
 *
 * \param t pointer to T - the tuple over which to "iterate" (elements will be
 * passed to the functor)
 * \param f F - the functor who shall be called with \a t and \a r as arguments
 * \param seq<Is...> - a struct that contains the sequence of numbers in its
 * template information
 */
template <typename T, typename F, int... Is> void for_each(T *t, F f, seq<Is...>) {
  auto l = {(f(&std::get<Is>(*t), Is), 0)...};
  (void)l; // suppress compiler warning about unused variable
}

/**
 * \brief Evaluate the functor for each element of the tuple. Not to be called
 * directly.
 *
 * This template also allows passing an additional argument to the functor.
 *
 * \param t pointer to T - the tuple over which to "iterate" (elements will be
 * passed to the functor)
 * \param r pointer to R - an argument that will be passed to the functor
 * \param f F - the functor who shall be called with \a t and \a r as arguments
 * \param seq<Is...> - a struct that contains the sequence of numbers in its
 * template information
 */
template <typename T, typename R, typename F, int... Is>
void for_each_vector(T *t, R *r, F f, seq<Is...>) {
  auto l = {(f(&std::get<Is>(*t), r, Is), 0)...};
  (void)l; // suppress compiler warning about unused variable
}

/**
 * \brief Evaluate the functor for each element of the tuple. Can be called
 * directly.
 *
 * This template also allows passing an additional argument to the functor.
 *
 * \param t pointer to tuple - the tuple over which to "iterate" (elements will be
 * passed to the functor)
 * \param r pointer to R - an argument that will be passed to the functor
 * \param f F - the functor who shall be called with \a t and \a r as arguments
 */
template <typename... Ts, typename R, typename F>
void for_each_in_tuple_vector(std::tuple<Ts...> *t, R *r, F f) {
  for_each_vector(t, r, f, gen_seq<sizeof...(Ts)>());
}

/**
 * \brief Evaluate the functor for each element of the tuple. Can be called
 * directly.
 *
 * \param t pointer to tuple - the tuple over which to "iterate" (elements will be
 * passed to the functor)
 * \param f F - the functor who shall be called with \a t and \a r as arguments
 */
template <typename... Ts, typename F>
void for_each_in_tuple(std::tuple<Ts...> *t, F f) {
  for_each(t, f, gen_seq<sizeof...(Ts)>());
}

} // namespace tuple_it

#endif // ITERATE_OVER_TUPLE_DEAMON_H
