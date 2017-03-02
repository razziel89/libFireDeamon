## \brief High level function that wraps the generation of a skin surface.
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
        raise ValueError("Lengths of coordinate list and radii list are not equal.")

    if shrink_factor<=0.0 or shrink_factor>=1.0:
        raise ValueError("Shrink factor must be between 0.0 and 1.0, excluding these borders.")

    if refinesteps<0:
        raise ValueError("The number of refinement steps must be a positive integer or 0.")
    
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

## \brief High level function that wraps the computation of the electrostatic potential via multithreaded C++ code.
def ElectrostaticPotentialPy(points, charges, coordinates, prog_report=True,cutoff=10000000.0):
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
        raise ValueError("Lengths of coordinate list and charges list are not equal.")

    charges_coordinates_vec=VectorDouble([cc for cc in _generate_three_one(coordinates,charges)]);
    points_vec=VectorDouble(list(iterchain.from_iterable(points)))

    potential_vec=VectorDouble()
    potential_vec.reserve(len(points))

    electrostatic_potential(prog_report, len(points), points_vec, charges_coordinates_vec, potential_vec,cutoff);
    
    potential=[p for p in potential_vec]

    del points_vec
    del charges_coordinates_vec
    del potential_vec
    
    return potential

## \brief Initialize data required to perform some computations on a grid
def InitializeGridCalculationOrbitalsPy(grid,basis,scale=1.0,normalize=True):
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
    scale: float, optional (default: 1.0)
        Divide each coordinate by this value (coordinate transformation).
    normalize: bool, optional (default: True)
        Whether or not to assume that the Gaussian functions that make up the
        primnitives are normalized or not.
    """
    import numpy as np
    vec_density_grid      = VectorDouble([p/scale for gpoint in grid for p in gpoint])
    density_indices       = np.array([ bc for b,bc in zip(basis,xrange(len(basis))) for prim in xrange(len(b[2])) ])
    vec_prim_centers      = VectorDouble([ p for b in basis for prim in b[2] for p in b[0] ])
    vec_prim_angular      = VectorInt([ a for b in basis for prim in b[2] for a in b[1] ])
    vec_prim_exponents    = VectorDouble([ prim[0] for b in basis for prim in b[2] ])
    vec_prim_coefficients = VectorDouble([ prim[1] for b in basis for prim in b[2] ])

    if normalize:
        normalize_gaussians(vec_prim_coefficients,vec_prim_exponents,vec_prim_angular)

    return vec_prim_centers, vec_prim_exponents, vec_prim_coefficients, vec_prim_angular, vec_density_grid, density_indices

## \brief Compute the electron density due to molecular orbitals
def ElectronDensityPy(coefficients_list,data,occupations=None,volume=1.0,prog_report=True,detailed_prog=False, cutoff=-1.0,correction=None):
    """
    Calculate the electron density due to some molecular orbitals on a grid.

    coefficients_list: list of lists of floats
        The coefficients of the molecular orbitals.
    data: what InitializeGridCalculationOrbitalsPy returned

    volume: float
        Scale the whole density by the inverse of this value.
    prog_report: bool
        Whether or not to give progress reports over MOs.
    detailed_prog:
        Whether or not to give progress reports while a MO
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

        tempdens =  np.array([d for d in density_vec])
        sumtemp  = np.sum(tempdens)
        if correction is not None:
            correction.append(1.0/sumtemp)
        if occupations is not None:
            density += occupations[rcount]*tempdens/sumtemp
        else:
            density += tempdens/sumtemp

        del density_vec
        del vec_mo_coefficients

        rcount+=1
        if prog_report:
            reportstring = "  %4d/%d"%(rcount,maxrcount)
            print ERASE_LINE+reportstring+CURSOR_UP_ONE

    if prog_report:
        print

    density*=1.0/volume

    density[density<1E-30]=0.0 #cut very small values away
    
    return density

## \brief High level wrapper to create an isosurface of arbitrary high discretization through volumetric data
def IsosurfacePy(data,origin,counts,delta,isovalue,points_inside,relative_precision=1.0e-05,mesh_criteria=[30,5,5]):
    """
    High level wrapper to create an isosurface of arbitrary high discretization
    through volumetric data. The data is given on an implicit regular grid in 3
    dimensions. One isosurface per element of points_inside is computed and
    overlaps are discarded. Using few points for points_inside greatly speeds
    up the computation.

    WARNING: if points_inside does not fit the data, the algorithm might not
             yield the actual iso surface.
    WARNING: the first mesh criterion (angular bound) is <30.0, the algorithm
             is not guaranteed to finish.
    HINT: if creating an iso-density-surface around a molecule, it is usually
          sufficient to pass the poisition of only one atom via points_inside.

    data: list of N floats
        A flat list of the volumetric data. The order of indices is that of
        dx-files, which is as follows:
            z - fast
            y - middle
            x - slow
    origin: list of 3 floats
        The origin of the 3 dimensional regular grid.
    counts: list of 3 int
        The number of points in each of the three directions of the grid.
        The product of these three values is the length of 'data'.
    delta: a 3x3 matrix (list of 3 lists with 3 elements each)
        The three vectors stored in this parameter form the vertex of the
        regular grid on which the data is defined. The matrix than can be
        built from these vectors must have any values unequal 0.0 solely
        on its main diagonal. This means that the three axes of the grid
        have to be aligned parallel to the three Cartesian axes.
    isovalue: float
        The isovalue at which to compute the isosurface.
    points_inside: an iterable of [float,float,float]
        Points that are expected to lie inside of (or at least very close to)
        the resulting isosurface. In the case of molecules, this can be the
        atoms' coordinates.
    relative_precision: float, optional (default: 1.0e-05)
        Precision value used to compute the isosurface. A lower value results
        in more highly discretized surfaces.
    mesh_criteria: a list of A,R,D, all floats. optional (default: [30.0,5.0,5.0]
        Explanations from: http://doc.cgal.org/latest/Surface_mesher/index.html
        A: float
            Angular bound for surface mesh generation. If <30, the algorithm is
            not guaranteed to finish. This is the lower bound in degrees for
            the angles during mesh generation.
        R: float
            Radius bound used during mesh generation. It is an upper bound on
            the radii of surface Delaunay balls. A surface Delaunay ball is a
            ball circumscribing a mesh facet and centered on the surface.
        D: float
            Distance bound used during surface mesh generation. It is an upper
            bound for the distance between the circumcenter of a mesh facet and
            the center of a surface Delaunay ball of this facet.
    """

    import numpy as np

    origin = np.array(origin)
    counts = np.array(counts)
    delta  = np.array(delta)

    if not np.allclose(delta, np.diag(np.diag(delta))):
        raise ValueError("ERROR: the voxel vectors must be parallel to the three cartesian axes.")

    delta_total = np.diag(delta)*counts

    data_vec   = VectorDouble([d for d in data]);
    extent_vec = VectorInt(list(counts))
    origin_vec = VectorDouble(list(origin))
    voxel_vec  = VectorDouble(list(np.diag(delta)))

    #points at all corners
    #points_inside=[origin+np.dot(delta_total,d) for d in np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]],dtype=float)]

    radii = []

    for point_inside in points_inside:
        radius = 0.0
        p = np.array(point_inside,dtype=float)
        for d in np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]],dtype=float):
            temp = np.linalg.norm(p-(origin+np.dot(delta_total,d)))
            if temp > radius:
                radius = temp
        radii.append(radius)

    radii_vec = VectorDouble(radii)
    points_inside_vec = VectorDouble([c for p in points_inside for c in p])

    mesh_criteria_vec = VectorDouble([float(e) for e in mesh_criteria])

    #for v in [origin_vec, voxel_vec, extent_vec,points_inside_vec, mesh_criteria_vec, radii_vec]:
    #    print [i for i in v]

    index_vec=VectorInt()
    length_vec=VectorInt()
    vertex_vec=VectorDouble()
    normal_vec=VectorDouble()

    make_isosurface(data_vec, origin_vec, voxel_vec, extent_vec,
            points_inside_vec, mesh_criteria_vec, radii_vec, relative_precision,
            isovalue,index_vec,vertex_vec,normal_vec,length_vec)

    nr_indices=length_vec.pop()
    nr_vertices=length_vec.pop()

    result=[(nr_indices,nr_vertices),
        [facet for facet in _generate_triples(index_vec,3*nr_indices)],
        [vertex for vertex in _generate_triples(vertex_vec,3*nr_vertices)],
        [n for n in _generate_triples(normal_vec,3*nr_vertices)]
    ]

    return result

## \brief Calculate the electron density due to some molecular orbitals on a grid.
def ElectrostaticPotentialOrbitalsPy(coefficients_list,Smat,occupations,data,prog_report=True):
    """
    Calculate the electron density due to some molecular orbitals on a grid.

    coefficients_list: list of lists of floats
        The coefficients of the molecular orbitals.
    Smat: list of lists of floats
        The overlap matrix between the primitive Gaussian functions
    occupations: a list of floats 
        The occupation number of the corresponding molecular orbital
    data: what InitializeGridCalculationOrbitalsPy returned

    prog_report: bool
        Whether or not to give progress reports over the grid.
    """
    import numpy as np

    vec_prim_centers, vec_prim_exponents, vec_prim_coefficients, vec_prim_angular, vec_potential_grid, potential_indices = data

    if not (vec_potential_grid.size())%3 == 0:
        raise ValueError("Number of values in vector for potential grid is not divisible by 3.")

    nr_gridpoints = int(vec_potential_grid.size()/3)

    nr_primitives = len(potential_indices)
    nr_electrons  = sum(occupations)
    np_o = np.array(occupations,dtype=float)
    np_c = np.array(coefficients_list,dtype=float).T
    #compute the first order density matrix P_mu_nu
    P = [
            [
                np.sum( np_o * np_c[potential_indices[mu]] * np_c[potential_indices[nu]] )
            for nu in xrange(nr_primitives)
            ]
        for mu in xrange(nr_primitives)
        ]

    potential_vec = VectorDouble()
    potential_vec.reserve(nr_gridpoints)

    P_vec = VectorDouble([p for pinner in P for p in pinner])
    #screen_vec = VectorInt([s for sinner in screen for s in sinner])

    S = [
            [
                Smat[potential_indices[mu]][potential_indices[nu]]
            for nu in xrange(nr_primitives)
            ]
        for mu in xrange(nr_primitives)
        ]
    S_vec = VectorDouble([s for sinner in S for s in sinner])

    electrostatic_potential_orbitals(prog_report, nr_gridpoints, vec_prim_centers, vec_prim_exponents, vec_prim_coefficients, vec_prim_angular, vec_potential_grid, P_vec, S_vec, potential_vec );

    potential = np.array([p for p in potential_vec])

    del potential_vec
    del P_vec

    return potential
