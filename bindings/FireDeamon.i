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
%}

%pythoncode %{
class LengthDisagreementError(Exception):
    pass

def _generate_coord_radii(coordinates,radii):
    for c,r in zip(coordinates,radii):
        yield c[0]
        yield c[1]
        yield c[2]
        yield r

def _generate_triples(l,length):
    for i in range(0,length,3):
        yield [l[i],l[i+1],l[i+2]]

def SkinSurfacePy(shrink_factor,coordinates,radii,refine=True,refinesteps=1):
    """
    High level function that wraps the generation of a skin surface.
    
    shrink_factor: shrink factor for the skin surface generation
    coordinates
    radii: a list containing all the radii 
    refine: whether or not to perform refinement (NOT YET IMPLEMENTED, ALWAYS True)
    refinesteps: refinement steps to perform (NOT YET IMPLEMENTED)
    """
    coord_radii_vec=VectorDouble([cr for cr in _generate_coord_radii(coordinates,radii)]);
    if len(coordinates)!=len(radii):
        raise LengthDisagreementError("Lengths of coordinate list and radii list are not equal.")
    
    nr_atoms=len(radii)
    
    vertex_vec=VectorDouble();
    index_vec=VectorInt();
    length_vec=VectorInt();
    
    make_skin_surface(shrink_factor,nr_atoms,coord_radii_vec,index_vec,vertex_vec,length_vec)
    nr_indices=length_vec.pop()
    nr_vertices=length_vec.pop()
    del length_vec
    del coord_radii_vec
    
    result=[(nr_indices,nr_vertices),[facet for facet in _generate_triples(index_vec,3*nr_indices)],[vertex for vertex in _generate_triples(vertex_vec,3*nr_vertices)]]
    
    del vertex_vec
    del index_vec
    
    return result
%}

%include "skin_surface_deamon.h"
