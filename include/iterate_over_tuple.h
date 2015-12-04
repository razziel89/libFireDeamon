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
#ifndef ITERATE_OVER_TUPLE_DEAMON_H
#define ITERATE_OVER_TUPLE_DEAMON_H

#include <tuple>

//taken from http://stackoverflow.com/questions/16387354/template-tuple-calling-a-function-on-each-element
namespace detail
{
    //generate a sequence of numbers
    template<int... Is>
    struct seq { };

    template<int N, int... Is>
    struct gen_seq : gen_seq<N - 1, N - 1, Is...> { };

    template<int... Is>
    struct gen_seq<0, Is...> : seq<Is...> { };

    //evaluate the functor for each element of the tuple
    //two parameters will be passed: the element of the tuple and where the result shall be stored (i.e. r)
    template<typename T, typename R, typename F, int... Is>
    void for_each_vector(T* t, R* r, F f, seq<Is...>)
    {
        auto l = { (f(&std::get<Is>(*t),r, Is), 0)... };
        (void) l; //suppress compiler warning about unused variable
    }

    //only the tuple elements will be passed over
    template<typename T, typename F, int... Is>
    void for_each(T* t, F f, seq<Is...>)
    {
        auto l = { (f(&std::get<Is>(*t), Is), 0)... };
        (void) l; //suppress compiler warning about unused variable
    }
}

template<typename... Ts, typename R, typename F>
void for_each_in_tuple_vector(std::tuple<Ts...>* t, R* r, F f)
{
    detail::for_each_vector(t, r, f, detail::gen_seq<sizeof...(Ts)>());
}

template<typename... Ts, typename F>
void for_each_in_tuple(std::tuple<Ts...>* t, F f)
{
    detail::for_each(t, f, detail::gen_seq<sizeof...(Ts)>());
}

#endif //ITERATE_OVER_TUPLE_DEAMON_H
