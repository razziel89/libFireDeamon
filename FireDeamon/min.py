## \package FireDeamon
#  \brief Python module for libFireDeamon
"""
This module includes Python wrapper functions to easily access the included
C++ functionality of libFireDeamon. It is simply called FireDeamon.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys
from itertools import chain as iterchain
from .cpp_min import *

def _generate_three_one(coordinates,radii):
    for c,r in zip(coordinates,radii):
        yield c[0]
        yield c[1]
        yield c[2]
        yield r

def _generate_triples(l,length):
    for i in range(0,length,3):
        yield [l[i],l[i+1],l[i+2]]

## \brief High level function that wraps the interpolation of arbitrary data on an irregular grid via multithreaded C++ code.
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
        raise ValueError("Lengths of coordinate list and vals list are not equal.")

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

    vals_vec=VectorDouble([float(v) for v in vals]);
    points_vec=VectorDouble([float(p) for p in iterchain.from_iterable(points)])
    coordinates_vec=VectorDouble([float(c) for c in iterchain.from_iterable(coordinates)])

    interpolation_vec=VectorDouble()
    interpolation_vec.reserve(len(points))
    
    generic_interpolation(prog_report, len(points), coordinates_vec, vals_vec, points_vec, interpolation_vec, interpolation_type, distance_exponent, distance_function, cutoff)

    interpolation=[i for i in interpolation_vec]

    del points_vec
    del coordinates_vec
    del vals_vec
    del interpolation_vec
    
    return interpolation 

## \brief deprecated
def NeighbourListPy(grid, nr_neighbours, cutoff, max_nr_neighbours=None, prog_report=True, cutoff_type='eukledian', sort_it=False):
    """
    Deprecated version of IrregularNeighbourListPy. Will be removed soon.
    """
    print("WARNING: use of NeighbourListPy is deprecated, use the new IrregularNeighbourListPy instead (same interface).",file=sys.stderr)
    return IrregularNeighbourListPy(grid, nr_neighbours, cutoff, max_nr_neighbours, prog_report, cutoff_type, sort_it)

## \brief Generate a list of neighbours of each point on an arbitrary grid
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
        print ("WARNING: not sorting with respect to distances but using max_nr_neighbours > nr_neighbours might throw away actual nearest neighbours.",file=sys.stderr)

    neighbours_vec = VectorInt()
    neighbours_vec.reserve(nr_gridpoints*(nr_neighbours+1));

    grid_vec = VectorDouble([float(component) for point in grid for component in point])

    make_neighbour_list_irregular(prog_report, nr_gridpoints, max_nr_neighbours, nr_neighbours, cutoff_int, grid_vec, cutoff_vec, neighbours_vec, sort_it);

    del grid_vec
    del cutoff_vec

    return neighbours_vec

## \brief Generate a list of neighbours of each point on a regular grid.
def RegularNeighbourListPy(nr_gridpoints_xyz, nr_neighbour_shells, prog_report=True, exclude_border=False):
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
    exclude_border: boolean, optional, default: False
        Whether or not points that do not have the maximum number of neighbours (i.e., those
        that lie close to or at the border) should be possible candidates for minima.
    """

    try:
        if not len(nr_gridpoints_xyz) == 3:
            raise ValueError("The variable nr_gridpoints_xyz must contain three values.")
    except TypeError as e:
        raise TypeError("The variable nr_gridpoints_xyz must be a list with three ints.",e)

    nr_gridpoints = 1
    for npts in nr_gridpoints_xyz:
        nr_gridpoints *= int(npts)
    if type(nr_gridpoints)!=int:
        raise TypeError("Too many points in grid, total number cannot be represented as integers.")

    nr_neighbours = (2*nr_neighbour_shells+1)**3 - 1

    neighbours_vec = VectorInt()
    neighbours_vec.reserve(nr_gridpoints*(nr_neighbours+1));

    make_neighbour_list_regular(prog_report, exclude_border, int(nr_gridpoints_xyz[0]), int(nr_gridpoints_xyz[1]), int(nr_gridpoints_xyz[2]), nr_neighbour_shells, neighbours_vec)

    return neighbours_vec

## \brief Given a neighbour list (as created by NeighbourListPy), find local minima
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
    values_vec = VectorDouble([float(v) for v in values])

    try:
        nr_degen = len(degeneration)
    except TypeError:
        degeneration = [degeneration]
        nr_degen = 1

    degeneration_vec = VectorDouble([float(d) for d in degeneration])

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
