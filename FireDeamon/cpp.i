// cpp.i
%module cpp

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

#ifndef FD_FULL_SUPPORT

%pythoncode %{
FULL_SUPPORT = False
%}

#else

%{
#include "FireDeamon/core/skin_surface_deamon.h"
#include "FireDeamon/core/isosurface.h"
#include "FireDeamon/core/electrostatic_potential_charges.h"
#include "FireDeamon/core/electrostatic_potential_orbitals.h"
#include "FireDeamon/core/electron_density.h"
#include "FireDeamon/core/orbital_overlap.h"
%}

%include "FireDeamon/core/skin_surface_deamon.h"
%include "FireDeamon/core/isosurface.h"
%include "FireDeamon/core/electrostatic_potential_charges.h"
%include "FireDeamon/core/electrostatic_potential_orbitals.h"
%include "FireDeamon/core/electron_density.h"
%include "FireDeamon/core/orbital_overlap.h"

%pythoncode %{
FULL_SUPPORT = True
%}

#endif

