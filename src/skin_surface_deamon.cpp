/***********
This file is part of libFireDeamon.

Copyright (C) 2015,2016 by Torsten Sachse

libFireDeamon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libFireDeamon is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with libFireDeamon.  If not, see <http://www.gnu.org/licenses/>.
***********/
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
#include <CGAL/Skin_surface_refinement_policy_3.h>
#include <list>
#include <vector>
#include <skin_surface_deamon.h>

template <class SkinSurface, class Polyhedron>
//extract data into 2 arrays and return length in 2 integers
//Since this is only to be used with polyhedrons without half-edges,
//that information will not be returned.
void get_polyhedron(SkinSurface &skin,
				   Polyhedron &p,
				   std::vector<int> &ivec,
                   std::vector<double> &dvec,
                   int *nr_vertices,
                   int *nr_facets)
{
    //define types
    typedef typename Polyhedron::Vertex_iterator                  Vertex_iterator;
    typedef typename Polyhedron::Facet_iterator                   Facet_iterator;
    typedef typename Polyhedron::Vertex_handle                    Vertex_handle;
    typedef typename Polyhedron::Halfedge_around_facet_circulator HFC;
    typedef typename Polyhedron::Traits::Point_3                  Point;


    CGAL::Skin_surface_refinement_policy_3<SkinSurface, Polyhedron> policy(skin);

    //get size of vectors
    *nr_vertices   = p.size_of_vertices ();
    *nr_facets     = p.size_of_facets   ();
//    int nr_halfedges  = p.size_of_halfedges()
    
    //resize vectors properly
    ivec.reserve(3 * *nr_facets  );
    dvec.reserve(3 * *nr_vertices);

    // add vertices to vector
    for (Vertex_iterator vit = p.vertices_begin(); vit != p.vertices_end(); ++vit) {
        Point vit_point = vit->point();
        dvec.push_back(vit_point.cartesian(0));
        dvec.push_back(vit_point.cartesian(1));
        dvec.push_back(vit_point.cartesian(2));
    }

    // add facets to vector
    CGAL::Inverse_index<Vertex_handle> index(p.vertices_begin(), p.vertices_end());
    for(Facet_iterator fi = p.facets_begin(); fi != p.facets_end(); ++fi) {
        HFC hc = fi->facet_begin();
        HFC hc_end = hc;
        std::size_t n = circulator_size( hc);
        do {
            Vertex_handle vh = (*hc).vertex();
            ivec.push_back(index[vh]);
        } while (++hc != hc_end);
    }
}

//declare type definitions
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Skin_surface_traits_3<K>                      Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::FT                                  FT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Weighted_point::Point                               Bare_point;
typedef CGAL::Polyhedron_3<K,
  CGAL::Skin_surface_polyhedral_items_3<Skin_surface_3> >   Polyhedron;

void make_skin_surface(double shrink_factor, int nr_atoms, std::vector<double> coord_radii_vec, std::vector<int> *ivec, std::vector<double> *dvec, std::vector<int> *length) {
    //declare variables for computation
    FT shrinkfactor = shrink_factor;
    std::list<Weighted_point> l;
    Polyhedron p;

    //read in centers of points from the input variables
    for (std::vector<double>::iterator it=coord_radii_vec.begin(); it!=coord_radii_vec.end();){
        double x,y,z,w;
        x=*it++;
        y=*it++;
        z=*it++;
        w=*it++;
        l.push_front(Weighted_point(Bare_point(x,y,z),w));
    }

    int *length_vert = (int*)malloc(sizeof(int));
    int *length_face = (int*)malloc(sizeof(int));
    length->reserve(2);
    
    //create the skin surface
    Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);
    CGAL::mesh_skin_surface_3(skin_surface, p);

    //subdivide, i.e., refine it
    CGAL::subdivide_skin_surface_mesh_3(skin_surface, p);
    
    //extract data from the generated skin_surface and the polyhedron
    //into 2 vectors and get their respective lengths
    
    //fill variables
    get_polyhedron(skin_surface, p, *ivec, *dvec, length_vert, length_face);
    length->push_back(*length_vert);
    length->push_back(*length_face);

    free(length_vert);
    free(length_face);

}

