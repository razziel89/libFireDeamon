/***********
This file is part of libFireDeamon.

Copyright (C) 2016 by Torsten Sachse

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
// Includes
// CGAL includes
#ifdef FDDEBUG
#define CGAL_DISABLE_ROUNDING_MATH_CHECK
#endif
#include "CGAL/Nef_polyhedron_3.h"
#include "CGAL/Polyhedron_3.h"
#include "CGAL/Polyhedron_incremental_builder_3.h"
#include "CGAL/Simple_cartesian.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
//
//---CGAL IO includes---
//#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
// use modified versions of the Image headers to remove dependencies on libGL.so and
// libGLU.so
#include <FireDeamon/core/CGAL/ImageIO_impl.h>
#include <FireDeamon/core/CGAL/Image_3_impl.h>
//#include <CGAL/ImageIO.h>
//
//---other default includes---
#include <deque>
#include <math.h>
#include <signal.h>
#include <stdexcept>
#include <vector>
//
//---deamon includes (interface)---
#include <FireDeamon/core/isosurface.h>

// triLinInterp

//---CGAL typedefs---
// Default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3; // c2t3
typedef Tr::Geom_traits GT;
typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;
typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kinex;
typedef CGAL::Exact_predicates_exact_constructions_kernel Kex;
typedef CGAL::Nef_polyhedron_3<Kex> Nef_Polyhedron;
typedef CGAL::Polyhedron_3<Kex> Ex_Polyhedron;
typedef CGAL::Polyhedron_3<Kinex> Polyhedron;

int signum = -1;
void iso_signal_callback_handler(int signal) {
  if (signum == -1) {
    fprintf(stderr,
            "\nCaught signal %d, interrupting at the next possible point\n",
            signal);
    signum = signal;
  } else {
    if (signal == signum) {
      fprintf(
          stderr, "\nCaught signal %d again, interrupting now (might fail)\n", signal);
      exit(signal);
    } else {
      fprintf(stderr,
              "\nCaught new signal %d, interrupting at the next possible point\n",
              signal);
      signum = signal;
    }
  }
}

/**
 * \brief A class that allows to copy a polyhedron declared on one kernel to a
 * polyhedron declared on another kernel.
 */
template <class Polyhedron_input, class Polyhedron_output>
struct Copy_polyhedron_to
    : public CGAL::Modifier_base<typename Polyhedron_output::HalfedgeDS> {
  Copy_polyhedron_to(const Polyhedron_input &in_poly) : in_poly(in_poly) {}

  void operator()(typename Polyhedron_output::HalfedgeDS &out_hds) {
    typedef typename Polyhedron_output::HalfedgeDS Output_HDS;

    CGAL::Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);

    typedef typename Polyhedron_input::Vertex_const_iterator Vertex_const_iterator;
    typedef typename Polyhedron_input::Facet_const_iterator Facet_const_iterator;
    typedef typename Polyhedron_input::Halfedge_around_facet_const_circulator HFCC;

    builder.begin_surface(in_poly.size_of_vertices(),
                          in_poly.size_of_facets(),
                          in_poly.size_of_halfedges());

    for (Vertex_const_iterator vi = in_poly.vertices_begin(),
                               end = in_poly.vertices_end();
         vi != end;
         ++vi) {
      typename Polyhedron_output::Point_3 p(::CGAL::to_double(vi->point().x()),
                                            ::CGAL::to_double(vi->point().y()),
                                            ::CGAL::to_double(vi->point().z()));
      builder.add_vertex(p);
    }

    typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
    Index index(in_poly.vertices_begin(), in_poly.vertices_end());

    for (Facet_const_iterator fi = in_poly.facets_begin(), end = in_poly.facets_end();
         fi != end;
         ++fi) {
      HFCC hc = fi->facet_begin();
      HFCC hc_end = hc;
      //     std::size_t n = circulator_size( hc);
      //     CGAL_assertion( n >= 3);
      builder.begin_facet();
      do {
        builder.add_vertex_to_facet(index[hc->vertex()]);
        ++hc;
      } while (hc != hc_end);
      builder.end_facet();
    }
    builder.end_surface();
  } // end operator()(..)
private:
  const Polyhedron_input &in_poly;
}; // end Copy_polyhedron_to<>

/**
 * \brief A wrapper function that allows to copy a polyhedron declared on one kernel to
 * a polyhedron declared on another kernel.
 *
 * This function wraps Copy_polyhedron_to for easy access. Be warned: not all
 * destination kernels can reproduce the same polyhedra that can be expressed using the
 * original kernel. Some might be fully incompatible.
 *
 * \param poly_a Poly_A - polyhedron on original kernel
 * \param poly_b Poly_B - polyhedron on destination kernel
 */
template <class Poly_A, class Poly_B>
void copy_to(const Poly_A &poly_a, Poly_B &poly_b) {
  Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
  poly_b.delegate(modifier);
}

/**
 * \brief A class for vector operations.
 */
struct Point3d;
typedef struct Point3d {

  /*! double - the point's x-coordinate */
  double x;
  /*! double - the point's y-coordinate */
  double y;
  /*! double - the point's z-coordinate */
  double z;
  /**
   * \brief Alternate constructor.
   *
   * When given a pointer to a double, take what thsi pointer points to
   * as the x-coordinate and the 2 values after that in memory as y- and
   * z-coordinates.
   *
   * \param p pointer to double - pointer to x-coordinate
   */
  Point3d(double *p) {
    x = *(p + 0);
    y = *(p + 1);
    z = *(p + 2);
  }
  //! \brief Default constructor.
  //!
  //! The point is initialized to the origin.
  Point3d() : x(0.0), y(0.0), z(0.0) {}
  /**
   * \brief Alternate constructor.
   *
   * \param _x double - x-coordinate
   * \param _y double - y-coordinate
   * \param _z double - z-coordinate
   */
  Point3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
  /**
   * \brief Copy constructor.
   *
   * \param p Point3d - point to copy
   */
  Point3d(const Point3d &p) : x(p.x), y(p.y), z(p.z) {}
  //! \brief subtract a vector
  struct Point3d operator-(struct Point3d p) {
    return Point3d(x - p.x, y - p.y, z - p.z);
  }
  //! \brief add a vector
  struct Point3d operator+(struct Point3d p) {
    return Point3d(x + p.x, y + p.y, z + p.z);
  }
  //! \brief add a vector directly
  struct Point3d &operator+=(const struct Point3d p) {
    x += p.x;
    y += p.y;
    z += p.z;
    return *this;
  }
  //! \brief subtract a vector directly
  struct Point3d &operator-=(const struct Point3d p) {
    x -= p.x;
    y -= p.y;
    z -= p.z;
    return *this;
  }
  //! \brief access the vector's 3 elements
  //! \return an element of the vector
  double operator[](int i) {
    if (i == 0) {
      return x;
    } else {
      if (i == 1) {
        return y;
      } else {
        if (i == 2) {
          return z;
        } else {
          throw std::logic_error("Index to point must be >=0 and <=2");
        }
      }
    }
  }
  //! \brief Compute the cross product of 2 vectors.
  //! \param p Point3d - vector with whom the cross product shall be computed
  struct Point3d operator*(struct Point3d p) {
    return Point3d(y * p.z - p.y * z, z * p.x - p.z * x, x * p.y - p.x * y);
  }
  //! \brief Scale a vector by a factor.
  //! \param d double - the scaling factor
  //! \return the scaled vector
  struct Point3d operator*(double d) {
    return Point3d(x * d, y * d, z * d);
  }
  //! \brief Scale a vector by the inverse of a factor.
  //! \param d double - the inverse of the scaling factor
  //! \return the scaled vector
  struct Point3d &operator/=(double d) {
    x /= d;
    y /= d;
    z /= d;
    return *this;
  }
  //! \brief Scale a vector by the inverse of a factor.
  //! \param i int - the inverse of the scaling factor
  //! \return the scaled vector
  struct Point3d &operator/=(unsigned int i) {
    double d = static_cast<double>(i);
    return *this /= d;
  }
  //! \brief Normalize this vector
  void normalize() {
    double norm = sqrt(x * x + y * y + z * z);
    x /= norm;
    y /= norm;
    z /= norm;
  }
} point;

std::ostream &operator<<(std::ostream &os, struct Point3d p) {
  os << "(" << p.x << " " << p.y << " " << p.z << ")";
  return os;
}

template <class Polyhedron>
// extract data into 2 arrays and return length in 2 integers
// Since this is only to be used with polyhedrons without half-edges,
// that information will not be returned.
void add_polyhedron(Polyhedron &p, std::vector<int> &ivec, std::vector<double> &dvec,
                    std::vector<double> &nvec, int *nr_vertices, int *nr_facets,
                    point origin) {
  // define types
  typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Halfedge_around_facet_circulator HFC;
  typedef typename Polyhedron::Traits::Point_3 Point;

  // double* dvec_first = dvec.end();
  // int*    ivec_first = ivec.end();
  int nr_points = *nr_vertices;

  // get size of vectors
  *nr_vertices += p.size_of_vertices();
  *nr_facets += p.size_of_facets();

  // reserve appropriate amount of memory
  ivec.reserve(3 * *nr_facets);
  dvec.reserve(3 * *nr_vertices);

  std::vector<std::pair<unsigned int, point>> normals;
  normals.reserve(p.size_of_vertices());

  // add vertices to vector
  for (Vertex_iterator vit = p.vertices_begin(); vit != p.vertices_end(); ++vit) {
    Point vit_point = vit->point();
    point p(vit_point.cartesian(0), vit_point.cartesian(1), vit_point.cartesian(2));
    p += origin;
    dvec.push_back(p.x);
    dvec.push_back(p.y);
    dvec.push_back(p.z);
    std::pair<unsigned int, point> temp_pair = std::make_pair(0, point());
    normals.push_back(temp_pair);
  }

  // add facets to vector
  CGAL::Inverse_index<Vertex_handle> index(p.vertices_begin(), p.vertices_end());
  for (Facet_iterator fi = p.facets_begin(); fi != p.facets_end(); ++fi) {
    HFC hc = fi->facet_begin();
    HFC hc_end = hc;
    std::size_t n = circulator_size(hc);
    if (n != 3) {
      char buff[1000];
      sprintf(buff, "Encountered facet with number of vertices unequal 3 but %lu.", n);
      throw std::logic_error(buff);
    }
    do {
      Vertex_handle vh = (*hc).vertex();
      ivec.push_back(nr_points + index[vh]);
    } while (++hc != hc_end);
    int i1 = *(ivec.end() - 3);
    int i2 = *(ivec.end() - 2);
    int i3 = *(ivec.end() - 1);
    point p1(&(dvec[3 * i1]));
    point p2(&(dvec[3 * i2]));
    point p3(&(dvec[3 * i3]));
    point normal(((p2 - p1) * (p3 - p2)));
    normal.normalize();
    normals[i1 - nr_points].first += 1;
    normals[i1 - nr_points].second += normal;
    normals[i2 - nr_points].first += 1;
    normals[i2 - nr_points].second += normal;
    normals[i3 - nr_points].first += 1;
    normals[i3 - nr_points].second += normal;
  }

  std::vector<std::pair<unsigned int, point>>::iterator normal_it;
  for (normal_it = normals.begin(); normal_it != normals.end(); ++normal_it) {
    if (normal_it->first == 0) {
      char buff[1000];
      sprintf(buff, "Found point that is in no triangle.");
      throw std::logic_error(buff);
    } else {
      normal_it->second /= normal_it->first;
      normal_it->second.normalize();
      nvec.push_back(normal_it->second.x);
      nvec.push_back(normal_it->second.y);
      nvec.push_back(normal_it->second.z);
    }
  }
}

bool intersect(std::deque<Nef_Polyhedron> &NefPolyDeq, Ex_Polyhedron &P,
               bool verbose = false, char verbose_string[] = NULL) {

  // if (not(P.is_closed())){
  //    throw std::logic_error("The computed polyhedron is not closed (which should be
  //    impossible).");
  //}
  bool intersection = false;
  Nef_Polyhedron NPin(P);
  if (verbose) {
    fprintf(stdout, "%c[2K\r", 27);
    fprintf(stdout, "%s - %-60s", verbose_string, "Converted to Nef-polyhedron");
    fflush(stdout);
  }
  for (std::deque<Nef_Polyhedron>::iterator it = NefPolyDeq.begin();
       !intersection && it != NefPolyDeq.end();
       ++it) {
    Nef_Polyhedron NPinter = *it * NPin;
    if (!NPinter.is_empty()) {
      intersection = true;
    }
  }
  if (verbose) {
    fprintf(stdout, "%c[2K\r", 27);
    fprintf(stdout, "%s - %-60s", verbose_string, "Checked for intersections");
    fflush(stdout);
  }
  if (!intersection) {
    NefPolyDeq.push_back(NPin);
  }
  return intersection;
}

template <typename T>
_image *get_image(std::vector<double> vec, std::vector<double> voxel, point center,
                  std::vector<int> extent) {

  _image *image = new _image();
  image->xdim = extent[0];
  image->ydim = extent[1];
  image->zdim = extent[2];
  image->vdim = 1;
  image->vx = voxel[0];
  image->vy = voxel[1];
  image->vz = voxel[2];
  image->cx = center[0];
  image->cy = center[1];
  image->cz = center[2];
  image->tx = 0.0;
  image->ty = 0.0;
  image->tz = 0.0;
  image->rx = 0.0;
  image->ry = 0.0;
  image->rz = 0.0;
  image->spm_offset = 0.0;
  image->spm_scale = 0.0;
  // Size of data type in bytes (4 for float and 8 for doulbe onl my machine)
  image->wdim = sizeof(T);
  image->imageFormat = NULL;
  image->vectMode = VM_SCALAR;
  image->wordKind = WK_FLOAT;
  image->sign = SGN_UNKNOWN;
  image->user = NULL;
  image->nuser = 0;
  image->fd = NULL;
  image->openMode = OM_CLOSE;
  image->endianness = END_LITTLE;
  image->dataMode = DM_BINARY;
  if (image->xdim * image->ydim * image->zdim != vec.size()) {
    throw std::invalid_argument("Given data vector does not have the same number of "
                                "elements as specified by the extent vector.");
  }
  T *data = (T *)malloc(sizeof(T) * image->xdim * image->ydim * image->zdim);
  T *data_it = data;
  std::vector<double>::iterator vec_it = vec.begin();
  int ex = extent[0];
  int ey = extent[1];
  int ez = extent[2];
  for (int xc = 0; xc < ex; ++xc) {
    for (int yc = 0; yc < ey; ++yc) {
      for (int zc = 0; zc < ez; ++zc, ++vec_it) {
        *(data_it + xc + yc * (ex) + zc * (ex * ey)) = static_cast<T>(*vec_it);
      }
    }
  }
  // for (; vec_it!=vec.end(); ++vec_it, ++data_it){
  //    *data_it = static_cast<T>(*vec_it);
  //}
  image->data = data;

  return image;
}

void make_isosurface(std::vector<double> data, std::vector<double> origin,
                     std::vector<double> voxel, std::vector<int> extent,
                     std::vector<double> points_inside,
                     std::vector<double> mesh_criteria, std::vector<double> radii,
                     double relative_precision, double isovalue, std::vector<int> *ivec,
                     std::vector<double> *dvec, std::vector<double> *nvec,
                     std::vector<int> *length) {
  // perform some sanity checks
  if (mesh_criteria.size() != 3) {
    throw std::invalid_argument("The mesh_criteria vector must have a length of 3.");
  }
  if (origin.size() != 3) {
    throw std::invalid_argument("The origin vector must have a length of 3.");
  }
  if (voxel.size() != 3) {
    throw std::invalid_argument("The voxel vector must have a length of 3.");
  }
  if (extent.size() != 3) {
    throw std::invalid_argument("The extent vector must have a length of 3.");
  }
  if ((points_inside.empty()) || (points_inside.size() % 3 != 0) ||
      (points_inside.size() / 3 != radii.size())) {
    throw std::invalid_argument("The points_inside vector does not have three times as "
                                "many elements as the radii vector.");
  }
  {
    double min = data[0];
    double max = min;
    for (std::vector<double>::iterator it = data.begin(); it != data.end(); ++it) {
      if (*it < min) {
        min = *it;
      }
      if (*it > max) {
        max = *it;
      }
    }
    if (isovalue < min || isovalue > max) {
      char buff[1000];
      sprintf(buff, "Isovalue %f not in data interval [%f,%f]", isovalue, min, max);
      throw std::invalid_argument(buff);
    }
  }
  point origin_point(&(origin[0]));
  point delta_point(extent[0] * voxel[0], extent[1] * voxel[1], extent[2] * voxel[2]);
  point center_point = origin_point + delta_point * 0.5;

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(
      mesh_criteria[0], mesh_criteria[1], mesh_criteria[2]);

  std::deque<Nef_Polyhedron> NefPolyDeq;
  // std::deque<Polyhedron>     PolyDeq;
  std::vector<double>::iterator radius, point_inside;

  int nr_vertices = 0;
  int nr_facets = 0;

  // std::cerr << "here" << std::endl << std::flush;
  bool verbose = (radii.size() > 1);
  if (verbose) {
    std::cout << "Started computation of isosurface starting from multiple points."
              << std::endl;
  }
  int current_radius = 0;
  int nr_radii = radii.size();
  char verbose_string[1000] = {'\0'};

  signal(SIGINT, iso_signal_callback_handler);

  radius = radii.begin();
  point_inside = points_inside.begin();
  for (; radius != radii.end() && signum == -1;
       ++radius, point_inside += 3, ++current_radius) {
    if (verbose) {
      sprintf(verbose_string,
              "Computing isosurface, internal point %d/%d",
              current_radius,
              nr_radii);
      fprintf(stdout, "%c[2K\r", 27);
      fprintf(stdout, "%s - %-50s", verbose_string, "Starting");
      fflush(stdout);
    }
    // std::cerr << "inner here 1" << std::endl << std::flush;
    _image *raw_image = get_image<double>(data, voxel, center_point, extent);
    // std::cerr << "inner here 1.5" << std::endl << std::flush;
    Tr tr;         // 3D-Delaunay triangulation
    C2t3 c2t3(tr); // 2D-complex in 3D-Delaunay triangulation
    // the 'function' is a 3D gray level image
    CGAL::Image_3 image3(raw_image);
    Gray_level_image image(image3, static_cast<float>(isovalue));
    // std::cerr << "inner here 2" << std::endl << std::flush;
    if (verbose) {
      fprintf(stdout, "%c[2K\r", 27);
      fprintf(stdout, "%s - %-50s", verbose_string, "Created grey level image");
      fflush(stdout);
    }

    // Carefully choosen bounding sphere: the center must be inside the
    // surface defined by 'image' and the radius must be high enough so that
    // the sphere actually bounds the whole image.
    GT::Point_3 bounding_sphere_center(
        *(point_inside + 0), *(point_inside + 1), *(point_inside + 2));
    GT::FT bounding_sphere_squared_radius = (*radius) * (*radius) * 2.0;
    GT::Sphere_3 bounding_sphere(bounding_sphere_center,
                                 bounding_sphere_squared_radius);
    // std::cerr << "inner here 3" << std::endl << std::flush;
    // definition of the surface with relative precision
    Surface_3 surface(image, bounding_sphere, relative_precision);
    if (verbose) {
      fprintf(stdout, "%c[2K\r", 27);
      fprintf(stdout, "%s - %-50s", verbose_string, "Defined enclosing surface");
      fflush(stdout);
    }
    // std::cerr << "inner here 4" << std::endl << std::flush;
    // meshing surface, with the "manifold without boundary" algorithm
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
    // std::cerr << "inner here 5" << std::endl << std::flush;
    if (verbose) {
      fprintf(stdout, "%c[2K\r", 27);
      fprintf(stdout, "%s - %-50s", verbose_string, "Meshed surface");
      fflush(stdout);
    }

    // extract data from mesh
    Polyhedron P;
    CGAL::output_surface_facets_to_polyhedron(c2t3, P);
    if (verbose) {
      fprintf(stdout, "%c[2K\r", 27);
      fprintf(stdout, "%s - %-50s", verbose_string, "Converted to inexact polyhedron");
      fflush(stdout);
    }
    // std::cerr << "Treating point number " << radius-radii.begin() << std::endl <<
    // std::flush;
    bool inter;
    if (radii.size() == 1) {
      inter = false;
    } else {
      Ex_Polyhedron Pex;
      copy_to<Polyhedron, Ex_Polyhedron>(P, Pex);
      if (verbose) {
        fprintf(stdout, "%c[2K\r", 27);
        fprintf(
            stdout, "%s - %-50s", verbose_string, "Transferred to exact polyhedron");
        fflush(stdout);
      }
      inter = intersect(NefPolyDeq, Pex, verbose, verbose_string);
    }
    if (verbose) {
      fprintf(stdout, "%c[2K\r", 27);
      fprintf(stdout,
              "%s - %-50s",
              verbose_string,
              "Checked for intersections with previous polyhedra");
      fflush(stdout);
    }
    if (!inter) {
      // PolyDeq.push_back(P);
      // std::cerr << "Not intersecting" << std::endl << std::flush;
      add_polyhedron(P, *ivec, *dvec, *nvec, &nr_vertices, &nr_facets, origin_point);
      if (verbose) {
        fprintf(stdout, "%c[2K\r", 27);
        fprintf(stdout,
                "%s - %-50s",
                verbose_string,
                "No intersection found, added polyhedron");
        fflush(stdout);
      }
    } else {
      if (verbose) {
        fprintf(stdout, "%c[2K\r", 27);
        fprintf(stdout, "%s - %-50s", verbose_string, "Intersection found, not adding");
        fflush(stdout);
      }
    }
    // free(raw_image->data);
    // delete raw_image;
  }
  // std::cerr << "here" << std::endl << std::flush;
  if (signum != -1) {
    fprintf(stderr, "Interrupting\n");
    exit(signum);
  }
  if (verbose) {
    fprintf(stdout, "%c[2K\r", 27);
    fprintf(stdout, "Computed all %d points.\n", nr_radii);
    fflush(stdout);
  }

  length->reserve(2);
  length->push_back(nr_vertices);
  length->push_back(nr_facets);

  // FILE* out = fopen("complete.off","w");
  // fprintf(out,"OFF\n%d %d %d\n",nr_vertices, nr_facets, 0);
  // for (std::vector<double>::iterator it=dvec->begin(); it!=dvec->end(); it+=3){
  //    fprintf(out,"%f %f %f\n",*(it+0),*(it+1),*(it+2));
  //}
  // for (std::vector<int>::iterator it=ivec->begin(); it!=ivec->end(); it+=3){
  //    fprintf(out,"3 %d %d %d\n",*(it+0),*(it+1),*(it+2));
  //}
  // fflush(out);
}
