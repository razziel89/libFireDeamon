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
//CGAL includes
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/Cartesian.h>
#include "CGAL/Simple_cartesian.h" 
#include "CGAL/Polyhedron_3.h" 
#include "CGAL/Nef_polyhedron_3.h" 
#include "CGAL/Polyhedron_incremental_builder_3.h" 
//
//---CGAL IO includes---
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/ImageIO.h>
//#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
//#include "CGAL/IO/Polyhedron_iostream.h" 
//
//---default IO includes---
#include <fstream>
#include <ostream>
//
//---other default includes---
#include <stdexcept>
#include <vector>
#include <deque>
#include <stdlib.h>
#include <math.h>
//
//---deamon includes (interface)---
#include <isosurface.h>

//---CGAL typedefs---
//typedef CGAL::Simple_cartesian<double> Kernel; 
//typedef CGAL::Polyhedron_3<Kernel> Polyhedron; 
//typedef CGAL::Nef_polyhedron_3<Kernel> Nef_Polyhedron; 
//typedef Polyhedron::HalfedgeDS HalfedgeDS; 
//typedef Nef_Polyhedron::Vector_3  Vector_3; 
//typedef Nef_Polyhedron::Aff_transformation_3  Aff_transformation_3; 
typedef CGAL::Surface_mesh_default_triangulation_3 Tr; // default triangulation for Surface_mesher
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;   // c2t3
typedef Tr::Geom_traits GT;
typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;
typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kinex;
typedef CGAL::Exact_predicates_exact_constructions_kernel Kex;
//typedef CGAL::Exact_integer  NT;
//typedef CGAL::Extended_homogeneous<NT>  Kernel;
//typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
typedef CGAL::Nef_polyhedron_3<Kex>     Nef_Polyhedron;
typedef CGAL::Polyhedron_3<Kex>         Ex_Polyhedron;
//typedef Ex_Polyhedron::HalfedgeDS       Ex_HalfedgeDS; 
typedef CGAL::Polyhedron_3<Kinex>       Polyhedron;

template <class Polyhedron_input,
          class Polyhedron_output>
struct Copy_polyhedron_to
        : public CGAL::Modifier_base<typename Polyhedron_output::HalfedgeDS>
{
    Copy_polyhedron_to(const Polyhedron_input& in_poly)
        : in_poly(in_poly) {}

    void operator()(typename Polyhedron_output::HalfedgeDS& out_hds)
    {
        typedef typename Polyhedron_output::HalfedgeDS Output_HDS;

        CGAL::Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);

        typedef typename Polyhedron_input::Vertex_const_iterator Vertex_const_iterator;
        typedef typename Polyhedron_input::Facet_const_iterator  Facet_const_iterator;
        typedef typename Polyhedron_input::Halfedge_around_facet_const_circulator HFCC;

        builder.begin_surface(in_poly.size_of_vertices(),
                              in_poly.size_of_facets(),
                              in_poly.size_of_halfedges());

        for(Vertex_const_iterator
            vi = in_poly.vertices_begin(), end = in_poly.vertices_end();
            vi != end ; ++vi)
        {
            typename Polyhedron_output::Point_3 p(::CGAL::to_double( vi->point().x()),
                                                  ::CGAL::to_double( vi->point().y()),
                                                  ::CGAL::to_double( vi->point().z()));
            builder.add_vertex(p);
        }

        typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
        Index index( in_poly.vertices_begin(), in_poly.vertices_end());

        for(Facet_const_iterator
            fi = in_poly.facets_begin(), end = in_poly.facets_end();
            fi != end; ++fi)
        {
            HFCC hc = fi->facet_begin();
            HFCC hc_end = hc;
            //     std::size_t n = circulator_size( hc);
            //     CGAL_assertion( n >= 3);
            builder.begin_facet ();
            do {
                builder.add_vertex_to_facet(index[hc->vertex()]);
                ++hc;
            } while( hc != hc_end);
            builder.end_facet();
        }
        builder.end_surface();
    } // end operator()(..)
private:
    const Polyhedron_input& in_poly;
}; // end Copy_polyhedron_to<>

template <class Poly_A, class Poly_B>
void copy_to(const Poly_A& poly_a, Poly_B& poly_b)
{
        Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
        poly_b.delegate(modifier);
}

template <class Polyhedron>
//extract data into 2 arrays and return length in 2 integers
//Since this is only to be used with polyhedrons without half-edges,
//that information will not be returned.
void get_polyhedron(Polyhedron &p,
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
    //typedef typename Polyhedron::Traits::Vector_3                 Vector;

    //get size of vectors
    *nr_vertices   = p.size_of_vertices ();
    *nr_facets     = p.size_of_facets   ();
//    int nr_halfedges  = p.size_of_halfedges()

    //reserve appropriate amount of memory
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
        //std::size_t n = circulator_size( hc);
        do {
            Vertex_handle vh = (*hc).vertex();
            ivec.push_back(index[vh]);
        } while (++hc != hc_end);
    }
}

bool intersect(std::deque<Nef_Polyhedron> & NefPolyDeq, Ex_Polyhedron & P){
    
    if (not(P.is_closed())){ 
        throw std::logic_error("The computed polyhedron is not closed (which should be impossible).");
    }
    bool intersection = false;
    Nef_Polyhedron NPin(P); 
    for (std::deque<Nef_Polyhedron>::iterator it=NefPolyDeq.begin(); !intersection && it!=NefPolyDeq.end(); ++it){
        Nef_Polyhedron NPinter = *it * NPin; 
        if (!NPinter.is_empty()){
            intersection = true;
        }
    }
    if (!intersection){
        NefPolyDeq.push_back(NPin);
    }
    return intersection;
} 

//template <class C2t3>
//bool get_surface_facets(std::ostream& os,
//                   const C2t3& c2t3){
//
//    using CGAL::Surface_mesher::number_of_facets_on_surface;
//    int options = CGAL::Surface_mesher::IO_ORIENT_SURFACE;
//    
//    typedef typename C2t3::Triangulation Tr;
//    typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
//    typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
//    typedef typename Tr::Facet Facet;
//    typedef typename Tr::Edge Edge;
//    typedef typename Tr::Vertex_handle Vertex_handle;
//    
//    // Header.
//    const Tr& tr = c2t3.triangulation();
//    
//    bool success = true;
//    
//    CGAL::Surface_mesher::Write_to_OFF_file<Tr>
//      off(os, (options & CGAL::Surface_mesher::IO_VERBOSE) != 0);
//    
//    success &= off.write_header(tr.number_of_vertices(),
//                    c2t3.number_of_facets());
//    
//    CGAL_assertion(c2t3.number_of_facets() == number_of_facets_on_surface(tr));
//    
//    // Finite vertices coordinates.
//    std::map<Vertex_handle, int> V;
//    int inum = 0;
//    for(Finite_vertices_iterator vit = tr.finite_vertices_begin();
//        vit != tr.finite_vertices_end();
//        ++vit)
//    {
//      V[vit] = inum++;
//      success &= off.write_vertex(vit);
//    }
//    
//    success &= off.begin_facets();
//    
//    Finite_facets_iterator fit = tr.finite_facets_begin();
//    std::set<Facet> oriented_set;
//    std::stack<Facet> stack;
//    
//    typename Tr::size_type number_of_facets = c2t3.number_of_facets();
//    
//    CGAL_assertion_code(typename Tr::size_type nb_facets = 0; )
//    
//      while (oriented_set.size() != number_of_facets)
//      {
//    while ( fit->first->is_facet_on_surface(fit->second) == false ||
//        oriented_set.find(*fit) != oriented_set.end() ||
//    
//        oriented_set.find(c2t3.opposite_facet(*fit)) !=
//        oriented_set.end() )
//    {
//      ++fit;
//    }
//    oriented_set.insert(*fit);
//    stack.push(*fit);
//    while(! stack.empty() )
//    {
//      Facet f = stack.top();
//      stack.pop();
//      for(int ih = 0 ; ih < 3 ; ++ih) {
//        const int i1  = tr.vertex_triple_index(f.second, tr. cw(ih));
//        const int i2  = tr.vertex_triple_index(f.second, tr.ccw(ih));
//    
//        const typename C2t3::Face_status face_status
//          = c2t3.face_status(Edge(f.first, i1, i2));
//        if(face_status == C2t3::REGULAR) {
//          Facet fn = c2t3.neighbor(f, ih);
//          if (oriented_set.find(fn) == oriented_set.end()) {
//        if(oriented_set.find(c2t3.opposite_facet(fn)) == oriented_set.end())
//        {
//          oriented_set.insert(fn);
//          stack.push(fn);
//        }
//        else {
//          success = false; // non-orientable
//        }
//          }
//        }
//        else if(face_status != C2t3::BOUNDARY) {
//          success = false; // non manifold, thus non-orientable
//        }
//      } // end "for each neighbor of f"
//    } // end "stack non empty"
//      } // end "oriented_set not full"
//    
//    for(typename std::set<Facet>::const_iterator fit =
//      oriented_set.begin();
//    fit != oriented_set.end();
//    ++fit)
//    {
//      const typename Tr::Cell_handle cell = fit->first;
//      const int& index = fit->second;
//      const int index1 = V[cell->vertex(tr.vertex_triple_index(index, 0))];
//      const int index2 = V[cell->vertex(tr.vertex_triple_index(index, 1))];
//      const int index3 = V[cell->vertex(tr.vertex_triple_index(index, 2))];
//      success &= off.write_facet(index1, index2, index3);
//      CGAL_assertion_code(++nb_facets);
//    }
//    
//    CGAL_assertion(nb_facets == number_of_facets);
//    
//    success &= off.write_footer();
//    return success;
//}

template <typename T>
_image* get_image(std::vector<double> vec, std::vector<double> center, std::vector<double> voxel, std::vector<int> extent){
    
    _image* image = new _image();
    image->xdim=extent[0]; image->ydim=extent[1]; image->zdim=extent[2];
    image->vdim=1;
    image->vx=voxel[0];    image->vy=voxel[1];    image->vz=voxel[2];
    image->cx=center[0];   image->cy=center[1];   image->cz=center[2];
    image->tx=0.0;         image->ty=0.0;         image->tz=0.0;
    image->rx=0.0;         image->ry=0.0;         image->rz=0.0;
    image->spm_offset=0.0; image->spm_scale=0.0;
    image->wdim=sizeof(T); //size of data type in bytes (4 for float and 8 for doulbe onl my machine)
    image->imageFormat = NULL;
    image->vectMode    = VM_SCALAR;
    image->wordKind    = WK_FLOAT;
    image->sign        = SGN_UNKNOWN;
    image->user        = NULL ;
    image->nuser       = 0;
    image->fd          = NULL;
    image->openMode    = OM_CLOSE;;
    image->endianness  = END_LITTLE;
    image->dataMode    = DM_BINARY;
    if (image->xdim * image->ydim * image->zdim != vec.size()){
        throw std::invalid_argument("Given data vector does not have the same number of elements as specified by the extent vector.");
    }
    T* data = (T*) malloc(sizeof(T) * image->xdim * image->ydim * image->zdim);
    T* data_it = data;
    std::vector<double>::iterator vec_it = vec.begin();
    for (; vec_it!=vec.end(); ++vec_it, ++data_it){
        *data_it = static_cast<T>(*vec_it);
    }
    image->data = data;

    return image;
}

//void print_image(_image* i){
//    std::cout << i->xdim << " " << i->ydim << " " << i->zdim << " " << i->vdim << std::endl;
//    std::cout << i->vx<< " " << i->vy<< " " << i->vz<< std::endl;
//    std::cout << i->tx<< " " << i->ty<< " " << i->tz<< std::endl;
//    std::cout << i->rx<< " " << i->ry<< " " << i->rz<< std::endl;
//    std::cout << i->cx<< " " << i->cy<< " " << i->cz<< std::endl;
//    std::cout << i->spm_offset<< " " << i->spm_scale<< std::endl;
//    std::cout << i->data<< std::endl;
//    std::cout << i->wdim<< std::endl;
//    std::cout << i->imageFormat<< std::endl;
//    std::cout << i->vectMode<< std::endl;
//    std::cout << i->wordKind<< std::endl;
//    std::cout << i->sign<< std::endl;
//    std::cout << i->user<< std::endl;
//    std::cout << i->nuser<< std::endl;
//    std::cout << i->fd<< std::endl;
//    std::cout << i->openMode<< std::endl;
//    std::cout << i->endianness<< std::endl;
//    std::cout << i->dataMode<< std::endl << std::flush;
//
//    //float* d = (float*) i->data;
//    //for (unsigned int j=1; j<= i->xdim * i->ydim * i->zdim; ++j, ++d){
//    //    std::cout << *d << " ";
//    //    if (j%10==0){
//    //        std::cout << std::endl;
//    //    }
//    //}
//}
//
//void print_floats(float f1, float f2){
//    char* c1 = reinterpret_cast<char*>(&f1);
//    char* c2 = reinterpret_cast<char*>(&f2);
//    char temp1[1000] = {'\0'};
//    char* tstart = temp1;
//    for (unsigned int i = 0; i < sizeof(float); ++i){
//        char C1 = c1[i];
//        for (int j = 7; j >= 0; --j){
//            int t = ((C1 >> j) & 1);
//            if (t==0){
//                *tstart = '1';
//            }else{
//                *tstart = '0';
//            }
//            ++tstart;
//        }
//    }
//    printf("%s\t",temp1);
//    char temp2[1000] = {'\0'};
//    tstart = temp2;
//    for (unsigned int i = 0; i < sizeof(float); ++i){
//        char C2 = c2[i];
//        for (int j = 7; j >= 0; --j){
//            int t = ((C2 >> j) & 1);
//            if (t==0){
//                *tstart = '1';
//            }else{
//                *tstart = '0';
//            }
//            ++tstart;
//        }
//    }
//    printf("%s\n",temp2);
//}

//typedef std::pair<std::vector<double>,std::vector<int>> mesh_data;
//
//mesh_data reduce_isosurfaces(std::deque<mesh_data> surfaces){
//    return surfaces[0];
//}

struct Point3d;

typedef struct Point3d {
    double x,y,z;
    Point3d(double* p){
        x = *(p+0);
        y = *(p+1);
        z = *(p+2);
    }
    Point3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z){}
    Point3d(const Point3d& p) : x(p.x), y(p.y), z(p.z){}
    struct Point3d operator-(struct Point3d p){
        return Point3d(x-p.x,y-p.y,z-p.z);
    }
    struct Point3d operator+(struct Point3d p){
        return Point3d(x+p.x,y+p.y,z+p.z);
    }
    struct Point3d& operator+=(const struct Point3d p){
        x += p.x;
        y += p.y;
        z += p.z;
        return *this;
    }
    //double operator*(struct Point3d p){
    //    return x*p.x + y*p.y + z*p.z;
    //}
    struct Point3d operator*(struct Point3d p){
        return Point3d(
                y*p.z - p.y*z,
                z*p.x - p.z*x,
                x*p.y - p.x*y
                );
    }
    struct Point3d& operator/=(double d){
        x /= d;
        x /= y;
        x /= z;
        return *this;
    }
    void normalize(){
        double norm = sqrt(x*x + y*y + z*z);
        x /= norm;
        y /= norm;
        z /= norm;
    }
} point;

std::ostream& operator<<(std::ostream& os, struct Point3d p){
    os << "(" << p.x << " " << p.y << " " << p.z << ")";
    return os;
}

void make_isosurface(std::vector<double> data, std::vector<double> center, std::vector<double> voxel, std::vector<int> extent,
        std::vector<double> points_inside, std::vector<double> mesh_criteria, std::vector<double> radii, double relative_precision,
        double isovalue, std::vector<int> *ivec, std::vector<double> *dvec, std::vector<double> *nvec, std::vector<int> *length){
    // perform some sanity checks
    if (mesh_criteria.size()!=3){
        throw std::invalid_argument("The mesh_criteria vector must have a length of 3.");
    }
    if (center.size()!=3){
        throw std::invalid_argument("The center vector must have a length of 3.");
    }
    if (voxel.size()!=3){
        throw std::invalid_argument("The voxel vector must have a length of 3.");
    }
    if (extent.size()!=3){
        throw std::invalid_argument("The extent vector must have a length of 3.");
    }
    if ( (points_inside.empty()) || (points_inside.size()%3!=0) || (points_inside.size()/3 != radii.size()) ){
        throw std::invalid_argument("The points_inside vector does not have three times as many elements as the radii vector.");
    }

    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(mesh_criteria[0],
                                                       mesh_criteria[1],
                                                       mesh_criteria[2]);

    //std::deque<mesh_data> surfaces;
    std::deque<Nef_Polyhedron> NefPolyDeq;
    std::deque<Polyhedron>     PolyDeq;
    std::vector<double>::iterator radius, point_inside;

    //Polyhedron all_surfaces;

    radius = radii.begin();
    point_inside = points_inside.begin();
    for (; radius!=radii.end(); ++radius, point_inside+=3){

        _image* raw_image = get_image<double>(data, center, voxel, extent);
        Tr tr;            // 3D-Delaunay triangulation
        C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
        // the 'function' is a 3D gray level image
        CGAL::Image_3 image3(raw_image);
        Gray_level_image image(image3, static_cast<float>(isovalue));

        // Carefully choosen bounding sphere: the center must be inside the
        // surface defined by 'image' and the radius must be high enough so that
        // the sphere actually bounds the whole image.
        GT::Point_3 bounding_sphere_center(*(point_inside+0), *(point_inside+1), *(point_inside+2));
        GT::FT bounding_sphere_squared_radius = (*radius)*(*radius)*2.0;
        GT::Sphere_3 bounding_sphere(bounding_sphere_center,
                                     bounding_sphere_squared_radius);
        // definition of the surface with relative precision
        Surface_3 surface(image, bounding_sphere, relative_precision);
        // meshing surface, with the "manifold without boundary" algorithm
        CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());

        //std::vector<double> double_dummy(100,0.0);
        //std::vector<int> int_dummy(100,0);

        //surfaces.push_back( std::make_pair(double_dummy,int_dummy) );

        //extract data from mesh
        Polyhedron P;
        CGAL::output_surface_facets_to_polyhedron(c2t3, P);
        Ex_Polyhedron Pex;
        copy_to<Polyhedron,Ex_Polyhedron>(P,Pex);
        //printf("Treating %u\t",(unsigned)(radius-radii.begin()));
        //fflush(stdout);
        bool inter = intersect(NefPolyDeq,Pex);
        //printf("Intersecting? %s\n",inter ? "yes" : "no");
        //fflush(stdout);
        if (!inter){
            PolyDeq.push_back(P);
            //char temp[100] = {'\0'};
            //sprintf(temp,"out_%u.off",(unsigned)(radius-radii.begin()));
            //std::ofstream out(temp);
            ////get_surface_facets(out, c2t3);
            //CGAL::output_surface_facets_to_off(out, c2t3);
            //std::cout << "File: " << temp << "\tFinal number of points: " << tr.number_of_vertices() << "\n";
            //fflush(stdout);
        }
    }

    Nef_Polyhedron Nef_result;
    for (std::deque<Nef_Polyhedron>::iterator nef_it=NefPolyDeq.begin(); nef_it!=NefPolyDeq.end(); ++nef_it){
        Nef_result += *nef_it;
    }

    Ex_Polyhedron ex_result;
    Nef_result.convert_to_Polyhedron(ex_result);

    Polyhedron result;
    copy_to<Ex_Polyhedron,Polyhedron>(ex_result,result);

    //std::vector<int> ivec;
    //std::vector<double> dvec;
    int nr_vertices;
    int nr_facets;

    get_polyhedron<Polyhedron>(result,
                   *ivec,
                   *dvec,
                   &nr_vertices,
                   &nr_facets);

    length->reserve(2);
    length->push_back(nr_vertices);
    length->push_back(nr_facets);

    std::vector<point> normals;
    normals.reserve(nr_facets);
    std::vector<std::deque<int>> my_triangles(nr_vertices,std::deque<int>());
    //my_triangles.reserve(nr_vertices);

    for (int triangle=0; triangle<nr_facets/3; ++triangle){
        int i1 = (*ivec)[3*triangle+0];
        int i2 = (*ivec)[3*triangle+1];
        int i3 = (*ivec)[3*triangle+2];
        //std::cout << triangle << ": " << i1 << "," << i2 << "," << i3 << std::endl << std::flush;
        std::cout << i1 << " " << i2 << " " << i3 << std::endl << std::flush;
        point p1(&( (*dvec)[3*i1] ));
        point p2(&( (*dvec)[3*i2] ));
        point p3(&( (*dvec)[3*i3] ));
        //std::cout << p1 << "\t" << p2 << "\t" << p3 << std::endl << std::flush;
        point temp( ((p2-p1)*(p3-p2)) );
        //std::cout << temp << std::endl << std::flush;
        temp.normalize();
        //std::cout << temp << std::endl << std::flush;
        normals.push_back( temp );
        (my_triangles[i1]).push_back(triangle);
        (my_triangles[i2]).push_back(triangle);
        (my_triangles[i3]).push_back(triangle);
        //std::cout << "pushed back" << std::endl << std::flush;
    }
    for(std::vector<std::deque<int>>::iterator vec_it=my_triangles.begin(); vec_it!=my_triangles.end(); ++vec_it){
        point p(0.0,0.0,0.0);
        //std::cout << "Adding " << vec_it-my_triangles.begin() << std::endl;
        for (std::deque<int>::iterator deq_it=vec_it->begin(); deq_it!=vec_it->end(); ++deq_it){
            p += normals[*deq_it];
            //std::cout << p << "\t" << normals[*deq_it] << std::endl;
        }
        p /= vec_it->size();
        nvec->push_back(p.x);
        nvec->push_back(p.y);
        nvec->push_back(p.z);
    }

    //FILE* out = fopen("complete.off","w");
    //fprintf(out,"OFF\n%d %d %d\n",nr_vertices, nr_facets, 0);
    //for (std::vector<double>::iterator it=dvec.begin(); it!=dvec.end(); it+=3){
    //    fprintf(out,"%f %f %f\n",*(it+0),*(it+1),*(it+2));
    //}
    //for (std::vector<int>::iterator it=ivec.begin(); it!=ivec.end(); it+=3){
    //    fprintf(out,"3 %d %d %d\n",*(it+0),*(it+1),*(it+2));
    //}
    //fflush(out);

}
