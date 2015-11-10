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
#include "electron_density.h"
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
    
    make_skin_surface(shrink_factor,coord_radii_vec,index_vec,vertex_vec,length_vec,refinesteps)
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

    electrostatic_potential(prog_report, len(points), points_vec, charges_coordinates_vec, potential_vec);
    
    potential=[p for p in potential_vec]

    del points_vec
    del charges_coordinates_vec
    del potential_vec
    
    return potential

def InterpolationPy(coordinates, vals, points, config=None, prog_report=True):
    """
    High level function that wraps the interpolation of arbitrary data
    on an irregular grid via multithreaded C++ code.
    
    coordinates: a list of 3-element elements containing the Cartesian coordinates
                 at which the given values are localized
    vals: a list of values
    points: a list of 3-element elements containing the Cartesian coordinates
            at which to interpolate
    config: a dictionary of configuration values. Keys are:
                method: the interpolation method ('nearest', 'distance' (for weighted inverse distance))
                for distance, also are needed:
                    exponent: exponent for norm
                    function: 2 equals Eukledian norm, 3 equals the three norm
    prog_report: whether or not to get progress reports during the computation
                 (since it can take long)
    """
    if len(coordinates)!=len(vals):
        raise LengthDisagreementError("Lengths of coordinate list and vals list are not equal.")

    if config is None: #default to narest neighbour interpolation
        interpolation_type=1
        distance_exponent=1.0
        distance_function=1
    else:
        if config['method']=='nearest':
            interpolation_type=1
            distance_exponent=1.0
            distance_function=1
        elif config['method']=='distance':
            interpolation_type=2
            distance_exponent=float(config['exponent'])
            distance_function=int(config['function'])
        else:
            raise WrongInterpolationTypeError("Supported interpolation types are: 'nearest' and 'distance'.")

    vals_vec=VectorDouble(vals);
    points_vec=VectorDouble(list(iterchain.from_iterable(points)))
    coordinates_vec=VectorDouble(list(iterchain.from_iterable(coordinates)))

    interpolation_vec=VectorDouble()
    interpolation_vec.reserve(len(points))
    
    generic_interpolation(prog_report, len(points), coordinates_vec, vals_vec, points_vec, interpolation_vec, interpolation_type, distance_exponent, distance_function)

    interpolation=[i for i in interpolation_vec]

    del points_vec
    del coordinates_vec
    del vals_vec
    del interpolation_vec
    
    return interpolation 

def InitializeElectronDensityPy(grid,basis,scale=1.0):
    """
    Create data structures suitable for efficiently computing
    the elctron density on an arbitrary grid. Call this first
    and then ElectronDensityPy(coefficients_list,data) where data
    is what this function returns.

    grid: list of [float,float,float]
        The Cartesian coordinates of the grid
    basis: a list of [A,L,Prim]
           with
           A: a list of 3 floats
                The center of the contracted Cartesian Gaussian function
           L: a list of 3 ints
                The polynomial exponents of the contracted Cartesian Gaussian
           Prim: a list of [alpha,pre]
                with
                alpha: float
                    The exponential factor of the primitive Gaussian function
                pre: float
                    The contraction coefficient of the primitive Gaussian function
    scale: float
        Divide each coordinate by this value (coordinate transformation).
    """
    import numpy as np
    vec_density_grid      = VectorDouble([p/scale for gpoint in grid for p in gpoint])
    density_indices       = np.array([ bc for b,bc in zip(basis,xrange(len(basis))) for prim in xrange(len(b[2])) ])
    vec_prim_centers      = VectorDouble([ p for b in basis for prim in b[2] for p in b[0] ])
    vec_prim_angular      = VectorDouble([ a for b in basis for prim in b[2] for a in b[1] ])
    vec_prim_exponents    = VectorDouble([ prim[0] for b in basis for prim in b[2] ])
    vec_prim_coefficients = VectorDouble([ prim[1] for b in basis for prim in b[2] ])

    return vec_prim_centers, vec_prim_exponents, vec_prim_coefficients, vec_prim_angular, vec_density_grid, density_indices

def ElectronDensityPy(coefficients_list,data,volume=1.0,prog_report=True,detailed_prog=False):

    import numpy as np

    if not prog_report and detailed_prog:
        detailed_prog = False

    if not detailed_prog:
        CURSOR_UP_ONE = '\x1b[1A'
        ERASE_LINE = '\x1b[2K'
    else:
        CURSOR_UP_ONE = ''
        ERASE_LINE = ''

    vec_prim_centers, vec_prim_exponents, vec_prim_coefficients, vec_prim_angular, vec_density_grid, density_indices = data

    if not (vec_density_grid.size())%3 == 0:
        raise ValueError("Number of values in vector for density grid is not divisible by 3.")

    nr_gridpoints = int(vec_density_grid.size()/3)

    density = np.zeros((nr_gridpoints,),dtype=float)

    #np_coeffs = np.array(coefficients_list)

    maxrcount = len(coefficients_list)

    if prog_report:
        print "  %4d/%d"%(0,maxrcount)+CURSOR_UP_ONE
        rcount = 0

    #for coefficients in np_coeffs:
    for coefficients in coefficients_list:

        vec_mo_coefficients = VectorDouble(list(np.array(coefficients)[density_indices]))

        density_vec = VectorDouble()
        density_vec.reserve(nr_gridpoints)

        electron_density(detailed_prog, nr_gridpoints, vec_prim_centers, vec_prim_exponents, vec_prim_coefficients, vec_prim_angular, vec_density_grid, vec_mo_coefficients, density_vec);

        density += np.array([d for d in density_vec])

        del density_vec
        del vec_mo_coefficients

        if prog_report:
            rcount+=1
            reportstring = "  %4d/%d"%(rcount,maxrcount)
            print ERASE_LINE+reportstring+CURSOR_UP_ONE

    if prog_report:
        print

    density/=volume

    density[density<1E-30]=0.0 #cut very small values away
    
    return density

%}

%include "skin_surface_deamon.h"
%include "parallel_generic.h"
%include "electrostatic_potential.h"
%include "irregular_grid_interpolation.h"
%include "electron_density.h"
