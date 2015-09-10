// FireDeamon.i
%module FireDeamon

%include "typemaps.i"
%include "std_vector.i"

namespace std {
    %template(VectorDouble) vector<double>;
    %template(VectorInt) vector<int>;
};

%{
#include "skin_surface_deamon.h"
#include "parallel_generic.h"
#include "electrostatic_potential.h"
#include "irregular_grid_interpolation.h"
%}

%pythoncode %{
class LengthDisagreementError(Exception):
    pass

class ShrinkFactorError(Exception):
    pass

class NumberRefinementsError(Exception):
    pass

class WrongInterpolationTypeError(Exception):
    pass

def _generate_three_one(coordinates,radii):
    for c,r in zip(coordinates,radii):
        yield c[0]
        yield c[1]
        yield c[2]
        yield r

def _generate_triples(l,length):
    for i in range(0,length,3):
        yield [l[i],l[i+1],l[i+2]]

def SkinSurfacePy(shrink_factor,coordinates,radii,refinesteps=1):
    """
    High level function that wraps the generation of a skin surface.
    
    shrink_factor: shrink factor for the skin surface generation
    coordinates: a list of cartesian coordinates declaring the centers of
                 the spheres
    radii: a list containing all the radii 
    refinesteps: refinement steps to perform. 0 will turn it off.
    """
    if len(coordinates)!=len(radii):
        raise LengthDisagreementError("Lengths of coordinate list and radii list are not equal.")

    if shrink_factor<=0.0 or shrink_factor>=1.0:
        raise ShrinkFactorError("Shrink factor must be between 0.0 and 1.0, excluding these borders.")

    if refinesteps<0:
        raise NumberRefinementsError("The number of refinement steps must be a positive integer or 0.")
    
    coord_radii_vec=VectorDouble([cr for cr in _generate_three_one(coordinates,radii)]);

    nr_atoms=len(radii)
    
    vertex_vec=VectorDouble();
    index_vec=VectorInt();
    length_vec=VectorInt();
    
    make_skin_surface(shrink_factor,nr_atoms,coord_radii_vec,index_vec,vertex_vec,length_vec,refinesteps)
    nr_indices=length_vec.pop()
    nr_vertices=length_vec.pop()
    del length_vec
    del coord_radii_vec
    
    result=[(nr_indices,nr_vertices),[facet for facet in _generate_triples(index_vec,3*nr_indices)],[vertex for vertex in _generate_triples(vertex_vec,3*nr_vertices)]]
    
    del vertex_vec
    del index_vec
    
    return result

from itertools import chain as iterchain

def ElectrostaticPotentialPy(points, charges, coordinates, prog_report=True):
    """
    High level function that wraps the computation of the electrostatic potential via
    multithreaded C++ code.
    
    points: a list of 3-element elements containing the Cartesian coordinates
            at which to compute the potential
    charges: a list of charges at the coordinates
    coordinates: a list of 3-element elements containing the Cartesian coordinates
                 at which the previously given charges are localized
    prog_report: whether or not to get progress reports during the computation
                 (since it can take long)
    """
    if len(charges)!=len(coordinates):
        raise LengthDisagreementError("Lengths of coordinate list and charges list are not equal.")

    charges_coordinates_vec=VectorDouble([cc for cc in _generate_three_one(coordinates,charges)]);
    points_vec=VectorDouble(list(iterchain.from_iterable(points)))

    potential_vec=VectorDouble()
    potential_vec.reserve(len(points))

    electrostatic_potential(prog_report, len(points), len(charges), points_vec, charges_coordinates_vec, potential_vec);
    
    potential=[p for p in potential_vec]

    del points_vec
    del charges_coordinates_vec
    del potential_vec
    
    return potential

def InterpolationPy(coordinates, vals, points, int_type='distance', prog_report=True, distance_exponent=1.0, distance_function=2):
    """
    High level function that wraps the interpolation of arbitrary data
    on an irregular grid via multithreaded C++ code.
    
    coordinates: a list of 3-element elements containing the Cartesian coordinates
                 at which the given values are localized
    vals: a list of values
    points: a list of 3-element elements containing the Cartesian coordinates
            at which to interpolate
    int_type: type of interpolation. Possible values are: 'nearest', 'distance' (for weighted inverse distance)
    prog_report: whether or not to get progress reports during the computation
                 (since it can take long)
    distance_exponent: for type 'distance', exponent for norm
    distance_function: for type 'distance', 2 equals Eukledian norm, 3 equals the three norm. Must be >0 and <4.
    """
    if len(coordinates)!=len(vals):
        raise LengthDisagreementError("Lengths of coordinate list and vals list are not equal.")

    if int_type=='nearest':
        interpolation_type=1
    elif int_type=='distance':
        interpolation_type=2
    else:
        raise WrongInterpolationTypeError("Supported interpolation types are: 'nearest' and 'distance'.")

    vals_vec=VectorDouble(vals);
    points_vec=VectorDouble(list(iterchain.from_iterable(points)))
    coordinates_vec=VectorDouble(list(iterchain.from_iterable(coordinates)))

    interpolation_vec=VectorDouble()
    interpolation_vec.reserve(len(points))
    
    generic_interpolation(prog_report, len(coordinates), len(vals), len(points), coordinates_vec, vals_vec, points_vec, interpolation_vec, interpolation_type, float(distance_exponent), distance_function)

    interpolation=[i for i in interpolation_vec]

    del points_vec
    del coordinates_vec
    del vals_vec
    del interpolation_vec
    
    return interpolation 
%}

%include "skin_surface_deamon.h"
%include "parallel_generic.h"
%include "electrostatic_potential.h"
%include "irregular_grid_interpolation.h"
