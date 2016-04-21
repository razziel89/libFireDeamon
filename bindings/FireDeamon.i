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
#include "electrostatic_potential.h"
#include "irregular_grid_interpolation.h"
#include "arbitrary_grid_local_minima.h"
#include "electron_density.h"
%}

%pythoncode %{

import sys

class LengthDisagreementError(Exception):
    pass

class ShrinkFactorError(Exception):
    pass

class NumberRefinementsError(Exception):
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
    normal_vec=VectorDouble();
    index_vec=VectorInt();
    length_vec=VectorInt();
    
    make_skin_surface(shrink_factor,coord_radii_vec,index_vec,vertex_vec,normal_vec,length_vec,refinesteps)
    nr_indices=length_vec.pop()
    nr_vertices=length_vec.pop()
    del length_vec
    del coord_radii_vec
    
    result=[(nr_indices,nr_vertices),
        [facet for facet in _generate_triples(index_vec,3*nr_indices)],
        [vertex for vertex in _generate_triples(vertex_vec,3*nr_vertices)],
        [n for n in _generate_triples(normal_vec,3*nr_vertices)]
    ]
    
    del vertex_vec
    del index_vec
    del normal_vec
    
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
        distance_exponent=1
        distance_function=1
        cutoff=-1.0
    else:
        if config['method']=='nearest':
            interpolation_type=1
            distance_exponent=1
            distance_function=1
            cutoff=-1.0
        elif config['method']=='distance':
            interpolation_type=2
            distance_exponent=int(config['exponent'])
            distance_function=int(config['function'])
            cutoff=float(config['cutoff'])
        else:
            raise TypeError("Supported interpolation types are: 'nearest' and 'distance'.")

    vals_vec=VectorDouble(vals);
    points_vec=VectorDouble(list(iterchain.from_iterable(points)))
    coordinates_vec=VectorDouble(list(iterchain.from_iterable(coordinates)))

    interpolation_vec=VectorDouble()
    interpolation_vec.reserve(len(points))
    
    generic_interpolation(prog_report, len(points), coordinates_vec, vals_vec, points_vec, interpolation_vec, interpolation_type, distance_exponent, distance_function, cutoff)

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

def ElectronDensityPy(coefficients_list,data,volume=1.0,prog_report=True,detailed_prog=False, cutoff=-1.0):
    """
    Calculate the electron density due to some molecular orbitals on a grid.

    coefficients_list: list of lists of floats
        The coefficients of the molecular orbitals.
    data: what InitializeElectronDensityPy returned

    volume: float
        Scale the whole density by the inverse of this value.
    prog_report: bool
        Whether or not ti give progress reports over MOs.
    detailed_prog:
        Whether or not ti give progress reports while a MO
        is being treated.
    cutoff: float in units of the grid
        No density will be computed if the difference between the
        gridpoint and the center of the basis function is larger
        than this value.
    """

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

    maxrcount = len(coefficients_list)

    if prog_report:
        print "  %4d/%d"%(0,maxrcount)+CURSOR_UP_ONE
        rcount = 0

    for coefficients in coefficients_list:

        vec_mo_coefficients = VectorDouble(list(np.array(coefficients)[density_indices]))

        density_vec = VectorDouble()
        density_vec.reserve(nr_gridpoints)

        electron_density(detailed_prog, nr_gridpoints, vec_prim_centers, vec_prim_exponents, vec_prim_coefficients, vec_prim_angular, vec_density_grid, vec_mo_coefficients, density_vec, cutoff);

        density += np.array([d for d in density_vec])

        del density_vec
        del vec_mo_coefficients

        if prog_report:
            rcount+=1
            reportstring = "  %4d/%d"%(rcount,maxrcount)
            print ERASE_LINE+reportstring+CURSOR_UP_ONE

    if prog_report:
        print

    density*=1.0/volume

    density[density<1E-30]=0.0 #cut very small values away
    
    return density

def NeighbourListPy(grid, nr_neighbours, cutoff, max_nr_neighbours=None, prog_report=True, cutoff_type='eukledian', sort_it=False):
    """
    Deprecated version of IrregularNeighbourListPy. Will be removed soon.
    """
    print >>sys.stderr,"WARNING: use of NeighbourListPy is deprecated, use the new IrregularNeighbourListPy instead (same interface)."
    return IrregularNeighbourListPy(grid, nr_neighbours, cutoff, max_nr_neighbours, prog_report, cutoff_type, sort_it)

def IrregularNeighbourListPy(grid, nr_neighbours, cutoff, max_nr_neighbours=None, prog_report=True, cutoff_type='eukledian', sort_it=False):
    """
    Generate a list of neighbours of each point on an arbitrary grid.

    grid: list of [float,float,float]
        The coordinates of each point in the grid.
    nr_neighbours: int
        How many neighbours shall be seeked per gridpoint.
    cutoff: float or [float,float,float] (depending on cutoff_type)
        Declare the cutoff distance for the given cutoff_type in units
        of the grid.
    max_nr_neighbours: int, optional, default: nr_neighbours
        The maximum number of neighbours to be searched per gridpoint.
        This cannot be smaller than nr_neighbours. If the given number
        of neighbours has been found within the given cutoff, no further
        neighbours are being searched. So you might not get the nearest
        ones if this value is too small. Greatly impacts performance.
    prog_report: boolean, optional, default: True
        Whether or not information about the progress of the calculation
        should be printed to stdout.
    cutoff_type: string, optional, default: eukledian
        define how to determine whether a gridpoint is to far away from
        another to be considered its neighbour. Possible values:
        eukledian:
            The distance is the absolute value of the difference vector.
            Requires cutoff to be one float.
        manhattan_single:
            The sum of the distances in x,y and z directions is the
            distance. Requires cutoff to be one float.
        manhattan_multiple:
            Treat each Cartesian direction independently. Requires
            cutoff to be [float,float,float].
    sort_it: boolean, optional, default: False
        Whether or not the neighbours found should be sorted with
        respect to the distance to the given point in increasing order.
    """

    nr_gridpoints = len(grid)

    cutoff_vec = VectorDouble()
    if cutoff_type.lower() == 'eukledian':
        cutoff_vec.reserve(1)
        try:
            cutoff_vec.push_back(cutoff)
        except TypeError as e:
            raise TypeError("The variable 'cutoff' is not of the correct type. Must be float.",e)
        cutoff_int = 1
    elif cutoff_type.lower() == 'manhattan_multiple':
        cutoff_vec.reserve(3)
        try:
            cutoff_vec.push_back(cutoff[0])
            cutoff_vec.push_back(cutoff[1])
            cutoff_vec.push_back(cutoff[2])
        except TypeError as e:
            raise TypeError("The variable 'cutoff' is not of the correct type. Must be a list of 3 floats.",e)
        cutoff_int = 2
    elif cutoff_type.lower() == 'manhattan_single':
        cutoff_vec.reserve(1)
        try:
            cutoff_vec.push_back(cutoff)
        except TypeError as e:
            raise TypeError("The variable 'cutoff' is not of the correct type. Must be float.",e)
        cutoff_int = 3
    else:
        raise TypeError("Wrong cutoff type given.")

    if max_nr_neighbours is None:
        max_nr_neighbours = nr_neighbours
    else:
        if max_nr_neighbours < nr_neighbours:
            raise ValueError("The value for max_nr_neighbours cannot be smaller than nr_neighbours.")

    if not sort_it and max_nr_neighbours > nr_neighbours:
        print >>sys.stderr, "WARNING: not sorting with respect to distances but using max_nr_neighbours > nr_neighbours might throw away actual nearest neighbours."

    neighbours_vec = VectorInt()
    neighbours_vec.reserve(nr_gridpoints*(nr_neighbours+1));

    grid_vec = VectorDouble([component for point in grid for component in point])

    make_neighbour_list_irregular(prog_report, nr_gridpoints, max_nr_neighbours, nr_neighbours, cutoff_int, grid_vec, cutoff_vec, neighbours_vec, sort_it);

    del grid_vec
    del cutoff_vec

    return neighbours_vec

def RegularNeighbourListPy(nr_gridpoints_xyz, nr_neighbour_shells, prog_report=True):
    """
    Generate a list of neighbours of each point on a regular grid.
    Returns a std::vector<int> with (2*nr_neighbour_shells+1)**3-1 elements per gridpoint indicating
    how many neighbours there are and which ones are the neighbours. If fewer than the maximum
    number of gridpoints was found (e.g. because the point is at a corner), -1's will be added.

    nr_gridpoints_xyz: [int, int, int]
        The number of points in each of the three directions of the regular 3D-grid.
    nr_neighbour_shells: int
        How many neighbour shells shall be treated (i.e., consider all those points to be
        neighbours who lie inside a cuboid that is spanned by 2*nr_neighbour_shells times the
        vectors that make up the grid. That cuboid is centered around each point.)
    prog_report: boolean, optional, default: True
        Whether or not information about the progress of the calculation
        should be printed to stdout.
    """

    try:
        if not len(nr_gridpoints_xyz) == 3:
            raise ValueError("The variable nr_gridpoints_xyz must contain three values.")
    except TypeError as e:
        raise TypeError("The variable nr_gridpoints_xyz must be a list with three ints.",e)

    for val in nr_gridpoints_xyz:
        if not isinstance(val,int):
            raise TypeError("nr_gridpoints_xyz must contain integers.",e)

    try:
        nr_gridpoints = reduce(int.__mul__,nr_gridpoints_xyz)
    except TypeError as e:
        raise TypeError("Too many points in grid, total number cannot be represented as integers.",e)

    nr_neighbours = (2*nr_neighbour_shells+1)**3 - 1

    neighbours_vec = VectorInt()
    neighbours_vec.reserve(nr_gridpoints*(nr_neighbours+1));

    make_neighbour_list_regular(prog_report, nr_gridpoints_xyz[0], nr_gridpoints_xyz[1], nr_gridpoints_xyz[2], nr_neighbour_shells, neighbours_vec)

    return neighbours_vec

def LocalMinimaPy(neighbour_list, values, degeneration, nr_neighbours, prog_report=False, upper_cutoff=None, lower_cutoff=None, sort_it=1, depths=None):
    """
    Given a neighbour list (as created by NeighbourListPy), find local minima. This is
    done by comparing the data at each point to that of its neighbours. The point is
    a local minimum if its associated value is at least 'degeneration' lower than that
    of all its neighbours.

    neighbour_list: std::vector<int> (or SWIG proxy)
        A list of neighbours. The format is: N, N1, N2, N3, ... NM
        and this repeats for every point. M is equal to nr_neighbours
        and N is the number of actual neighbours that have been found
        for the respective point.
    values: list of floats
        The volumetric data in which the local minima shall be found.
    degeneration: float
        As mentioned in the above description. Can be positive or negative.
    prog_report: boolean, optional, default: False
        Whether or not information about the progress of the calculation
        should be printed to stdout.
    upper_cutoff: float, optional, default: do not use
        Do not consider points as possible minima whose associated values
        are at least this large.
    lower_cutoff: float, optional, default: do not use
        Do not consider points as possible minima whose associated values
        are at most this large.
    sort_it: int, optional, default: 1
        If 0, do not sort the minima by depths and do not return the depths.
        If 1, sort the minima with respect to the difference between the value
        at the minimum and the average of all surrounding points. If depths
        is not None, also append the estimated depths.
        If 2, sort the minima with respect to the difference between the value
        at the minimum and the minimum value of all surrounding points. If depths
        is not None, also append the estimated depths.
    depths: object that has an 'append' method
        if not None, append to the list (or other object) the estimated
        depths of the minima according to the value of sort_it. If it does
        not have this method, do not append the depths.
    """

    minima_vec = VectorInt()
    values_vec = VectorDouble(values)

    try:
        nr_degen = len(degeneration)
    except TypeError:
        degeneration = [degeneration]
        nr_degen = 1

    degeneration_vec = VectorDouble(degeneration)

    nr_values = len(values)

    use_upper = (upper_cutoff is not None)
    use_lower = (lower_cutoff is not None)

    if not use_upper:
        upper_cutoff = 0.0
    if not use_lower:
        lower_cutoff = 0.0

    depths_vec = None
    if depths is not None:
        if hasattr(depths,'append'):
            depths_vec = VectorDouble()

    if neighbour_list.size() != values_vec.size()*(nr_neighbours+1):
        raise ValueError("I expect %d ints in neighbour_list for each float in values, but I found %.2f."
                %( nr_neighbours+1, 1.0*neighbour_list.size()/values_vec.size() )
                )

    local_minima_from_neighbour_list(prog_report, nr_neighbours, nr_values, neighbour_list, values_vec, minima_vec, degeneration_vec, use_upper, use_lower, upper_cutoff, lower_cutoff, sort_it, depths_vec)

    del values_vec

    minima = [e for e in minima_vec]

    if depths_vec is not None:
        for d in depths_vec:
            depths.append(d)

    del depths_vec 

    return minima

%}

%include "skin_surface_deamon.h"
%include "electrostatic_potential.h"
%include "irregular_grid_interpolation.h"
%include "arbitrary_grid_local_minima.h"
%include "electron_density.h"
