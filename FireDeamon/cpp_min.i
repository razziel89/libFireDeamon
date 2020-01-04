// cpp_min.i
%module cpp_min

%include "typemaps.i"
%include "std_vector.i"
%include "std_string.i"

namespace std {
    %template(VectorDouble) vector<double>;
    %template(VectorInt) vector<int>;
};

%{
#include "FireDeamon/core/irregular_grid_interpolation.h"
#include "FireDeamon/core/arbitrary_grid_local_minima.h"
#include "FireDeamon/core/set_procname.h"
%}

%include "FireDeamon/core/irregular_grid_interpolation.h"
%include "FireDeamon/core/arbitrary_grid_local_minima.h"
%include "FireDeamon/core/set_procname.h"
