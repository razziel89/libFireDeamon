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
#include "FireDeamon/skin_surface_deamon.h"
#include "FireDeamon/isosurface.h"
#include "FireDeamon/electrostatic_potential_charges.h"
#include "FireDeamon/electrostatic_potential_orbitals.h"
#include "FireDeamon/irregular_grid_interpolation.h"
#include "FireDeamon/arbitrary_grid_local_minima.h"
#include "FireDeamon/electron_density.h"
#include "FireDeamon/orbital_overlap.h"
#include "FireDeamon/set_procname.h"
%}

%include "FireDeamon/skin_surface_deamon.h"
%include "FireDeamon/isosurface.h"
%include "FireDeamon/electrostatic_potential_charges.h"
%include "FireDeamon/electrostatic_potential_orbitals.h"
%include "FireDeamon/irregular_grid_interpolation.h"
%include "FireDeamon/arbitrary_grid_local_minima.h"
%include "FireDeamon/electron_density.h"
%include "FireDeamon/orbital_overlap.h"
%include "FireDeamon/set_procname.h"
