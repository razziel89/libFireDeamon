"""
This module includes Python wrapper functions to easily access the included
C++ functionality of libFireDeamon. It is simply called FireDeamon.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys
import os
from itertools import chain as iterchain
from .cpp import *

if os.environ.get("FD_REQ_FULL_SUPPORT", "0") == "1" and not FULL_SUPPORT:
    raise ImportError("Full support requested but unavailable")


def _generate_three_one(coordinates, radii):
    """Yield single elements from 2 iterables, 3 of the 1st one and 1 of the 2nd one

    Args:
        coordinates: (iterable) - after three elements of this iterable, one elements of
            the other one will be yielded, and then three of this one, and so on
        radii: (iterable) - after three elements of the other iterable, one elements of
            this one one will be yielded, and then three of the other one, and so on
    """
    for c, r in zip(coordinates, radii):
        yield c[0]
        yield c[1]
        yield c[2]
        yield r


def _generate_triples(l, length):
    """Convert an iterable intro triplets

    Args:
        l: (iterable) - iterable to be converted
        length: (int) - number of elements of the iterable
    """
    for i in range(0, length, 3):
        yield [l[i], l[i + 1], l[i + 2]]


def InterpolationPy(coordinates, vals, points, config=None, prog_report=True):
    """Interpolatee arbitrary data on an irregular grid via multithreaded C++ code.

    Args:
        coordinates: a list of 3-element elements containing the Cartesian coordinates
            at which the given values are localized

        vals: the list of given values

        points: a list of 3-element elements containing the Cartesian coordinates at
            which to interpolate

        config: a dictionary of configuration values, keys are

            method: the interpolation method 'nearest' for nearest neighbours,
            'distance' for weighted inverse distance

            exponent: exponent for norm (int) (only for 'distance')

            function: 2 equals Eukledian norm, 3 equals the three norm (int)  (only for
            'distance')

            cutoff: points farther away than this in the specified norm are ignored, set
            to a negative value to include all points (float) (only for 'distance')

        prog_report: whether or not to get progress reports during the computation
            (since it can take long) (only for 'distance')
    """
    if len(coordinates) != len(vals):
        raise ValueError("Lengths of coordinate list and vals list are not equal.")

    if config is None:  # default to narest neighbour interpolation
        interpolation_type = 1
        distance_exponent = 1
        distance_function = 1
        cutoff = -1.0
    else:
        if config["method"] == "nearest":
            interpolation_type = 1
            distance_exponent = 1
            distance_function = 1
            cutoff = -1.0
        elif config["method"] == "distance":
            interpolation_type = 2
            distance_exponent = int(config["exponent"])
            distance_function = int(config["function"])
            cutoff = float(config["cutoff"])
        else:
            raise TypeError(
                "Supported interpolation types are: 'nearest' and 'distance'."
            )

    # Convert all input lists to flattened vectors
    vals_vec = VectorDouble([float(v) for v in vals])
    points_vec = VectorDouble([float(p) for p in iterchain.from_iterable(points)])
    coordinates_vec = VectorDouble(
        [float(c) for c in iterchain.from_iterable(coordinates)]
    )

    interpolation_vec = VectorDouble()
    interpolation_vec.reserve(len(points))

    generic_interpolation(
        prog_report,
        len(points),
        coordinates_vec,
        vals_vec,
        points_vec,
        interpolation_vec,
        interpolation_type,
        distance_exponent,
        distance_function,
        cutoff,
    )

    # Convert result to a Python list
    interpolation = [i for i in interpolation_vec]

    del points_vec
    del coordinates_vec
    del vals_vec
    del interpolation_vec

    return interpolation


def IrregularNeighbourListPy(
    grid,
    nr_neighbours,
    cutoff,
    max_nr_neighbours=None,
    prog_report=True,
    cutoff_type="eukledian",
    sort_it=False,
):
    """ Generate a list of neighbours of each point on an arbitrary grid.

    This is a multi-step procedure. For each point, first, all other points are located
    that lie within the specified cutoff. Let N be that number. However, if N would be
    larger than max_nr_neighbours, only the first max_nr_neighbours many points will be
    retained and the other N-max_nr_neighbours many points will be silently ignored.
    Next, the first nr_neighbours many points out of those N points are found and
    retained. If sort_it is True, the N points are first sorted by distance before
    returning nr_neighbours many. This process is performed for eery grid point. If
    N<max_nr_neighbours and with a large cutoff, many points will be checked for their
    distance to the reference point, which can negatively impact performance.

    Args:
        grid: (list of [float,float,float]) - The coordinates of each point in the grid.
        nr_neighbours: (int) - How many neighbours shall be seeked per gridpoint.
        cutoff: (float or [float,float,float]) - (depending on cutoff_type) Declare the
            cutoff distance for the given cutoff_type in units of the grid.
        max_nr_neighbours: (int, optional, default: nr_neighbours) - The maximum number
            of neighbours to be searched per gridpoint.  This cannot be smaller than
            nr_neighbours. If the given number of neighbours has been found within the
            given cutoff, no further neighbours are being searched. So you might not get
            the nearest ones if this value is too small. Greatly impacts performance.
        prog_report: (boolean, optional, default: True) - Whether or not information
            about the progress of the calculation should be printed to stdout.
        cutoff_type: (string, optional, default: "eukledian")

            define how to determine whether a gridpoint is too far away from another to
            be considered its neighbour. Possible values:

            eukledian: The distance is the absolute value of the difference vector.
            Requires cutoff to be one float.

            manhattan_single: The sum of the distances in x,y and z directions is the
            distance. Requires cutoff to be one float.

            manhattan_multiple: Treat each Cartesian direction independently. Requires
            cutoff to be [float,float,float].

        sort_it: (boolean, optional, default: False) - Whether or not the neighbours
            found should be sorted with respect to the distance to the given point in
            increasing order.
    """

    nr_gridpoints = len(grid)

    cutoff_vec = VectorDouble()
    if cutoff_type.lower() == "eukledian":
        cutoff_vec.reserve(1)
        try:
            cutoff_vec.push_back(cutoff)
        except TypeError as e:
            raise TypeError(
                "The variable 'cutoff' is not of the correct type. Must be float.", e
            )
        cutoff_int = 1
    elif cutoff_type.lower() == "manhattan_multiple":
        cutoff_vec.reserve(3)
        try:
            cutoff_vec.push_back(cutoff[0])
            cutoff_vec.push_back(cutoff[1])
            cutoff_vec.push_back(cutoff[2])
        except TypeError as e:
            raise TypeError(
                "The variable 'cutoff' is not of the correct type. "
                "Must be a list of 3 floats.",
                e,
            )
        cutoff_int = 2
    elif cutoff_type.lower() == "manhattan_single":
        cutoff_vec.reserve(1)
        try:
            cutoff_vec.push_back(cutoff)
        except TypeError as e:
            raise TypeError(
                "The variable 'cutoff' is not of the correct type. Must be float.", e
            )
        cutoff_int = 3
    else:
        raise TypeError("Wrong cutoff type given.")

    if max_nr_neighbours is None:
        max_nr_neighbours = nr_neighbours
    else:
        if max_nr_neighbours < nr_neighbours:
            raise ValueError(
                "The value for max_nr_neighbours cannot be smaller than nr_neighbours."
            )

    if not sort_it and max_nr_neighbours > nr_neighbours:
        print(
            "WARNING: not sorting with respect to distances but using max_nr_neighbours"
            " > nr_neighbours might throw away actual nearest neighbours.",
            file=sys.stderr,
        )

    neighbours_vec = VectorInt()
    neighbours_vec.reserve(nr_gridpoints * (nr_neighbours + 1))

    grid_vec = VectorDouble([float(component) for point in grid for component in point])

    make_neighbour_list_irregular(
        prog_report,
        nr_gridpoints,
        max_nr_neighbours,
        nr_neighbours,
        cutoff_int,
        grid_vec,
        cutoff_vec,
        neighbours_vec,
        sort_it,
    )

    del grid_vec
    del cutoff_vec

    return neighbours_vec


def RegularNeighbourListPy(
    nr_gridpoints_xyz, nr_neighbour_shells, prog_report=True, exclude_border=False
):
    """Generate a list of neighbours of each point on a regular grid.

    Returns a std::vector<int> with (2*nr_neighbour_shells+1)**3-1 elements per
    gridpoint indicating how many neighbours there are and which ones are the
    neighbours. If fewer than the maximum number of gridpoints was found (e.g. because
    the point is at a corner), -1's will be added.

    Args:
        nr_gridpoints_xyz: ([int, int, int]) - The number of points in each of the three
            directions of the regular 3D-grid.
        nr_neighbour_shells: (int) - How many neighbour shells shall be treated (i.e.,
            consider all those points to be neighbours who lie inside a cuboid that is
            spanned by 2*nr_neighbour_shells times the vectors that make up the grid.
            That cuboid is centered around each point.)
        prog_report: (boolean, optional, default: True) - Whether or not information
            about the progress of the calculation should be printed to stdout.
        exclude_border: (boolean, optional, default: False) - Whether or not points that
            do not have the maximum number of neighbours (i.e., those that lie close to
            or at the border, depending on nr_neighbour_shells) should be considered to
            have neighbours. This is important when finding minima with LocalMinimaPy.
            When setting exclude_border to True, those excluded points are not
            considered possible minima.
    """

    try:
        if not len(nr_gridpoints_xyz) == 3:
            raise ValueError(
                "The variable nr_gridpoints_xyz must contain three values."
            )
    except TypeError as e:
        raise TypeError(
            "The variable nr_gridpoints_xyz must be a list with three ints.", e
        )

    nr_gridpoints = 1
    for npts in nr_gridpoints_xyz:
        nr_gridpoints *= int(npts)
    if type(nr_gridpoints) != int:
        raise TypeError(
            "Too many points in grid, total number cannot be represented as integers."
        )

    nr_neighbours = (2 * nr_neighbour_shells + 1) ** 3 - 1

    neighbours_vec = VectorInt()
    neighbours_vec.reserve(nr_gridpoints * (nr_neighbours + 1))

    make_neighbour_list_regular(
        prog_report,
        exclude_border,
        int(nr_gridpoints_xyz[0]),
        int(nr_gridpoints_xyz[1]),
        int(nr_gridpoints_xyz[2]),
        nr_neighbour_shells,
        neighbours_vec,
    )

    return neighbours_vec


def LocalMinimaPy(
    neighbour_list,
    values,
    degeneration,
    nr_neighbours,
    prog_report=False,
    upper_cutoff=None,
    lower_cutoff=None,
    sort_it=1,
    depths=None,
):
    """Search local minima.

    Given a neighbour list (as created by RegularNeighbourListPy or
    IrregularNeighbourListPy), find local minima. This is done by comparing the data at
    each point to that of its neighbours. The point is a local minimum if its associated
    value is at least 'degeneration' lower than that of all its neighbours.

    Args:
        neighbour_list: (std::vector<int> (or SWIG proxy)) - A list of neighbours. The
            format is: N, N1, N2, N3, ... NM and this repeats for every point. M is
            equal to nr_neighbours and N is the number of actual neighbours that have
            been found for the respective point.
        values: (list of floats) - The volumetric data in which the local minima shall
            be found.
        degeneration: (float) - As mentioned in the above description. Can be positive
            or negative.
        prog_report: (boolean, optional, default: False) - Whether or not information
            about the progress of the calculation should be printed to stdout.
        upper_cutoff: (float, optional, default: do not use) - Do not consider points as
            possible minima whose associated values are at least this large.
        lower_cutoff: (float, optional, default: do not use) - Do not consider points as
            possible minima whose associated values are at most this large.
        sort_it: (int, optional, default: 1)

            If 0, do not sort the minima by depths and do not return the depths.

            If 1, sort the minima with respect to the difference between the value at
            the minimum and the average of all surrounding points. If depths is not
            None, also append the estimated depths.

            If 2, sort the minima with respect to the difference between the value at
            the minimum and the minimum value of all surrounding points. If depths is
            not None, also append the estimated depths.

        depths: (object that has an 'append' method) - if not None, append to the list
            (or other object) the estimated depths of the minima according to the value
            of sort_it. If it does not have this method, do not append the depths.
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

    use_upper = upper_cutoff is not None
    use_lower = lower_cutoff is not None

    if not use_upper:
        upper_cutoff = 0.0
    if not use_lower:
        lower_cutoff = 0.0

    depths_vec = None
    if depths is not None:
        if hasattr(depths, "append"):
            depths_vec = VectorDouble()

    if neighbour_list.size() != values_vec.size() * (nr_neighbours + 1):
        raise ValueError(
            "I expect %d ints in neighbour_list for each float in values, but I found %.2f."
            % (nr_neighbours + 1, 1.0 * neighbour_list.size() / values_vec.size())
        )

    local_minima_from_neighbour_list(
        prog_report,
        nr_neighbours,
        nr_values,
        neighbour_list,
        values_vec,
        minima_vec,
        degeneration_vec,
        use_upper,
        use_lower,
        upper_cutoff,
        lower_cutoff,
        sort_it,
        depths_vec,
    )

    del values_vec

    minima = [e for e in minima_vec]

    if depths_vec is not None:
        for d in depths_vec:
            depths.append(d)

    del depths_vec

    return minima


# Define functions available with full support
if FULL_SUPPORT:

    def SkinSurfacePy(shrink_factor, coordinates, radii, refinesteps=1):
        """ High level function that wraps the generation of a skin surface.
        
        Args:
            shrink_factor: (float) - shrink factor for the skin surface generation
            coordinates: a list of cartesian coordinates declaring the centers of the
                spheres
            radii: a list containing all the radii 
            refinesteps: (int) - refinement steps to perform. 0 will turn it off.
        """
        if len(coordinates) != len(radii):
            raise ValueError("Lengths of coordinate list and radii list are not equal.")

        if shrink_factor <= 0.0 or shrink_factor >= 1.0:
            raise ValueError(
                "Shrink factor must be between 0.0 and 1.0, excluding these borders."
            )

        if refinesteps < 0:
            raise ValueError(
                "The number of refinement steps must be a positive integer or 0."
            )

        coord_radii_vec = VectorDouble(
            [float(cr) for cr in _generate_three_one(coordinates, radii)]
        )

        nr_atoms = len(radii)

        vertex_vec = VectorDouble()
        normal_vec = VectorDouble()
        index_vec = VectorInt()
        length_vec = VectorInt()

        make_skin_surface(
            shrink_factor,
            coord_radii_vec,
            index_vec,
            vertex_vec,
            normal_vec,
            length_vec,
            refinesteps,
        )
        nr_indices = length_vec.pop()
        nr_vertices = length_vec.pop()
        del length_vec
        del coord_radii_vec

        result = [
            (nr_indices, nr_vertices),
            [facet for facet in _generate_triples(index_vec, 3 * nr_indices)],
            [vertex for vertex in _generate_triples(vertex_vec, 3 * nr_vertices)],
            [n for n in _generate_triples(normal_vec, 3 * nr_vertices)],
        ]

        del vertex_vec
        del index_vec
        del normal_vec

        return result

    def ElectrostaticPotentialPy(
        points, charges, coordinates, prog_report=True, cutoff=10000000.0
    ):
        """Compute the electrostatic potential.

        High level function that wraps the computation of the electrostatic potential
        via multithreaded C++ code.
        
        Args:
            points: a list of 3-element elements containing the Cartesian coordinates at
                which to compute the potential
            charges: a list of charges at the coordinates
            coordinates: a list of 3-element elements containing the Cartesian
                coordinates at which the previously given charges are localized
            prog_report: whether or not to get progress reports during the computation
                (since it can take long)
        """
        if len(charges) != len(coordinates):
            raise ValueError(
                "Lengths of coordinate list and charges list are not equal."
            )

        charges_coordinates_vec = VectorDouble(
            [float(cc) for cc in _generate_three_one(coordinates, charges)]
        )
        points_vec = VectorDouble([float(p) for p in iterchain.from_iterable(points)])

        potential_vec = VectorDouble()
        potential_vec.reserve(len(points))

        electrostatic_potential(
            prog_report,
            len(points),
            points_vec,
            charges_coordinates_vec,
            potential_vec,
            cutoff,
        )

        potential = [p for p in potential_vec]

        del points_vec
        del charges_coordinates_vec
        del potential_vec

        return potential

    def InitializeGridCalculationOrbitalsPy(grid, basis, scale=1.0, normalize=True):
        """Initialize a grid calculation based on orbitals.

        Create data structures suitable for efficiently computing the elctron density on
        an arbitrary grid. Call this first and then
        ElectronDensityPy(coefficients_list,data) where data is what this function
        returns.
    
        Args:
            grid: (list of [float,float,float]) - The Cartesian coordinates of the grid
            basis: a list of [A,L,Prim] with

                   A: a list of 3 floats The center of the contracted Cartesian Gaussian
                   function

                   L: a list of 3 ints The polynomial exponents of the contracted
                   Cartesian Gaussian

                   Prim: a list of [alpha,pre] with

                        alpha: float The exponential factor of the primitive Gaussian
                        function

                        pre: float The contraction coefficient of the primitive Gaussian
                        function

            scale: (float, optional, default: 1.0) - Divide each coordinate by this
                value (coordinate transformation).
            normalize: (bool, optional, default: True) - Whether or not to assume that
                the Gaussian functions that make up the primnitives are normalized or
                not.
        """
        import numpy as np

        vec_density_grid = VectorDouble(
            [float(1.0 * p / scale) for gpoint in grid for p in gpoint]
        )
        density_indices = np.array(
            [bc for b, bc in zip(basis, range(len(basis))) for prim in range(len(b[2]))]
        )
        vec_prim_centers = VectorDouble(
            [float(p) for b in basis for prim in b[2] for p in b[0]]
        )
        vec_prim_angular = VectorInt(
            [int(a) for b in basis for prim in b[2] for a in b[1]]
        )
        vec_prim_exponents = VectorDouble(
            [float(prim[0]) for b in basis for prim in b[2]]
        )
        vec_prim_coefficients = VectorDouble(
            [float(prim[1]) for b in basis for prim in b[2]]
        )

        if normalize:
            normalize_gaussians(
                vec_prim_coefficients, vec_prim_exponents, vec_prim_angular
            )

        return (
            vec_prim_centers,
            vec_prim_exponents,
            vec_prim_coefficients,
            vec_prim_angular,
            vec_density_grid,
            density_indices,
        )

    def ElectronDensityPy(
        coefficients_list,
        data,
        occupations=None,
        volume=1.0,
        prog_report=True,
        detailed_prog=False,
        cutoff=-1.0,
        correction=None,
    ):
        """ Calculate the electron density due to some molecular orbitals on a grid.

        Args:
            coefficients_list: list of lists of floats The coefficients of the molecular
                orbitals.
            data: what InitializeGridCalculationOrbitalsPy returned
            volume: (float) - Scale the whole density by the inverse of this value.
            prog_report: (bool) - Whether or not to give progress reports over MOs.
            detailed_prog: (bool) - Whether or not to give progress reports while a MO
                is being treated.
            cutoff: (float in units of the grid) - No density will be computed if the
                difference between the gridpoint and the center of the basis function is
                larger than this value.
        """
        import numpy as np

        if not prog_report and detailed_prog:
            detailed_prog = False

        if not detailed_prog:
            CURSOR_UP_ONE = "\x1b[1A"
            ERASE_LINE = "\x1b[2K"
        else:
            CURSOR_UP_ONE = ""
            ERASE_LINE = ""

        (
            vec_prim_centers,
            vec_prim_exponents,
            vec_prim_coefficients,
            vec_prim_angular,
            vec_density_grid,
            density_indices,
        ) = data

        if not (vec_density_grid.size()) % 3 == 0:
            raise ValueError(
                "Number of values in vector for density grid is not divisible by 3."
            )

        nr_gridpoints = int(vec_density_grid.size() / 3)

        density = np.zeros((nr_gridpoints,), dtype=float)

        maxrcount = len(coefficients_list)

        if prog_report:
            print("  %4d/%d" % (0, maxrcount) + CURSOR_UP_ONE)
        rcount = 0

        for coefficients in coefficients_list:

            vec_mo_coefficients = VectorDouble(
                [float(moc) for moc in np.array(coefficients)[density_indices]]
            )

            density_vec = VectorDouble()
            density_vec.reserve(nr_gridpoints)

            electron_density(
                detailed_prog,
                nr_gridpoints,
                vec_prim_centers,
                vec_prim_exponents,
                vec_prim_coefficients,
                vec_prim_angular,
                vec_density_grid,
                vec_mo_coefficients,
                density_vec,
                cutoff,
            )

            tempdens = np.array([d for d in density_vec])
            sumtemp = np.sum(tempdens)
            if correction is not None:
                correction.append(1.0 / sumtemp)
            if occupations is not None:
                density += occupations[rcount] * tempdens / sumtemp
            else:
                density += tempdens / sumtemp

            del density_vec
            del vec_mo_coefficients

            rcount += 1
            if prog_report:
                reportstring = "  %4d/%d" % (rcount, maxrcount)
                print(ERASE_LINE + reportstring + CURSOR_UP_ONE)

        if prog_report:
            print()

        density *= 1.0 / volume

        density[density < 1e-30] = 0.0  # cut very small values away

        return density

    def IsosurfacePy(
        data,
        origin,
        counts,
        delta,
        isovalue,
        points_inside,
        relative_precision=1.0e-05,
        mesh_criteria=[30, 5, 5],
    ):
        """Generate an isosurface.

        High level wrapper to create an isosurface of arbitrary high discretization
        through volumetric data. The data is given on an implicit regular grid in 3
        dimensions. One isosurface per element of points_inside is computed and overlaps
        are discarded. Using few points for points_inside greatly speeds up the
        computation.
    
        Beware the following:

        - If points_inside does not fit the data, the algorithm might not
          yield the actual iso surface.
        - If the first mesh criterion (angular bound) is <30.0, the algorithm is not be
          guaranteed to finish.
        - If creating an iso-density-surface around a molecule, it is usually sufficient
          to pass the poisition of only one atom via points_inside.
    
        Args:
            data:
                list of N floats

                A flat list of the volumetric data. The order of indices is that of
                dx-files, which is as follows: z: fast, y: middle, x: slow

            origin: (list of 3 floats) - The origin of the 3 dimensional regular grid.
            counts: (list of 3 int) - The number of points in each of the three
                directions of the grid.  The product of these three values is the length
                of 'data'.
            delta: (a 3x3 matrix (list of 3 lists with 3 elements each)) - The three
                vectors stored in this parameter form the vertex of the regular grid on
                which the data is defined. The matrix than can be built from these
                vectors must have any values unequal 0.0 solely on its main diagonal.
                This means that the three axes of the grid have to be aligned parallel
                to the three Cartesian axes.
            isovalue: (float) - The isovalue at which to compute the isosurface.
            points_inside: (an iterable of [float,float,float]) - Points that are
                expected to lie inside of (or at least very close to) the resulting
                isosurface. In the case of molecules, this can be the atoms' coordinates.
            relative_precision: (float, optional, default: 1.0e-05) - Precision value
                used to compute the isosurface. A lower value results in more highly
                discretized surfaces.
            mesh_criteria:
                a list of A,R,D, all floats, optional ,default: [30.0,5.0,5.0]

                Explanations from: http://doc.cgal.org/latest/Surface_mesher/index.html

                A: Angular bound for surface mesh generation. If <30, the
                algorithm is not guaranteed to finish. This is the lower bound in
                degrees for the angles during mesh generation.

                R: Radius bound used during mesh generation. It is an upper bound on the
                radii of surface Delaunay balls. A surface Delaunay ball is a ball
                circumscribing a mesh facet and centered on the surface.

                D: Distance bound used during surface mesh generation. It is an upper
                bound for the distance between the circumcenter of a mesh facet and the
                center of a surface Delaunay ball of this facet.
        """

        import numpy as np

        origin = np.array(origin)
        counts = np.array(counts)
        delta = np.array(delta)

        if not np.allclose(delta, np.diag(np.diag(delta))):
            raise ValueError(
                "ERROR: the voxel vectors must be parallel to the three cartesian axes."
            )

        delta_total = np.diag(delta) * counts

        data_vec = VectorDouble([float(d) for d in data])
        extent_vec = VectorInt([int(c) for c in counts])
        origin_vec = VectorDouble([float(o) for o in origin])
        voxel_vec = VectorDouble([float(d) for d in np.diag(delta)])

        # points at all corners
        # points_inside=[origin+np.dot(delta_total,d) for d in np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]],dtype=float)]

        radii = []

        for point_inside in points_inside:
            radius = 0.0
            p = np.array(point_inside, dtype=float)
            for d in np.array(
                [
                    [0, 0, 0],
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1],
                    [1, 1, 0],
                    [1, 0, 1],
                    [0, 1, 1],
                    [1, 1, 1],
                ],
                dtype=float,
            ):
                temp = np.linalg.norm(p - (origin + np.dot(delta_total, d)))
                if temp > radius:
                    radius = temp
            radii.append(radius)

        radii_vec = VectorDouble([float(r) for r in radii])
        points_inside_vec = VectorDouble([float(c) for p in points_inside for c in p])

        mesh_criteria_vec = VectorDouble([float(e) for e in mesh_criteria])

        # for v in [origin_vec, voxel_vec, extent_vec,points_inside_vec, mesh_criteria_vec, radii_vec]:
        #    print([i for i in v])

        index_vec = VectorInt()
        length_vec = VectorInt()
        vertex_vec = VectorDouble()
        normal_vec = VectorDouble()

        make_isosurface(
            data_vec,
            origin_vec,
            voxel_vec,
            extent_vec,
            points_inside_vec,
            mesh_criteria_vec,
            radii_vec,
            relative_precision,
            isovalue,
            index_vec,
            vertex_vec,
            normal_vec,
            length_vec,
        )

        nr_indices = length_vec.pop()
        nr_vertices = length_vec.pop()

        result = [
            (nr_indices, nr_vertices),
            [facet for facet in _generate_triples(index_vec, 3 * nr_indices)],
            [vertex for vertex in _generate_triples(vertex_vec, 3 * nr_vertices)],
            [n for n in _generate_triples(normal_vec, 3 * nr_vertices)],
        ]

        return result

    def ElectrostaticPotentialOrbitalsPy(
        coefficients_list, Smat, occupations, data, prog_report=True
    ):
        """ Calculate the electron density due to some molecular orbitals on a grid.
    
        Args:
            coefficients_list: (list of lists of floats) - The coefficients of the
                molecular orbitals.
            Smat: (list of lists of floats) - The overlap matrix between the primitive
                Gaussian functions
            occupations: (a list of floats) - The occupation number of the corresponding
                molecular orbital
            data: what InitializeGridCalculationOrbitalsPy returned
            prog_report: (bool) - Whether or not to give progress reports over the grid.
        """
        import numpy as np

        (
            vec_prim_centers,
            vec_prim_exponents,
            vec_prim_coefficients,
            vec_prim_angular,
            vec_potential_grid,
            potential_indices,
        ) = data

        if not (vec_potential_grid.size()) % 3 == 0:
            raise ValueError(
                "Number of values in vector for potential grid is not divisible by 3."
            )

        nr_gridpoints = int(vec_potential_grid.size() / 3)

        nr_primitives = len(potential_indices)
        nr_electrons = sum(occupations)
        np_o = np.array(occupations, dtype=float)
        np_c = np.array(coefficients_list, dtype=float).T
        # compute the first order density matrix P_mu_nu
        primlist = list(range(nr_primitives))
        P = [
            [
                np.sum(np_o * np_c[potential_indices[mu]] * np_c[potential_indices[nu]])
                for nu in primlist
            ]
            for mu in primlist
        ]

        potential_vec = VectorDouble()
        potential_vec.reserve(nr_gridpoints)

        P_vec = VectorDouble([float(p) for pinner in P for p in pinner])
        # screen_vec = VectorInt([s for sinner in screen for s in sinner])

        S = [
            [Smat[potential_indices[mu]][potential_indices[nu]] for nu in primlist]
            for mu in primlist
        ]
        S_vec = VectorDouble([float(s) for sinner in S for s in sinner])

        electrostatic_potential_orbitals(
            prog_report,
            nr_gridpoints,
            vec_prim_centers,
            vec_prim_exponents,
            vec_prim_coefficients,
            vec_prim_angular,
            vec_potential_grid,
            P_vec,
            S_vec,
            potential_vec,
        )

        potential = np.array([p for p in potential_vec])

        del potential_vec
        del P_vec

        return potential
