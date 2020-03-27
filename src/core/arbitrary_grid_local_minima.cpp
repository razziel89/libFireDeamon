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
#include <FireDeamon/core/arbitrary_grid_local_minima.h>
#include <FireDeamon/core/parallel_generic.h>
#include <algorithm>
#include <limits>
#include <math.h>
#include <pthread.h>
#include <stdexcept>
#include <time.h>
#include <tuple>
#include <utility>
#include <vector>

bool sort_by_first(std::pair<double, int> p1, std::pair<double, int> p2) {
  return p1.first < p2.first;
}

double get_average(int nr_values, std::vector<int>::iterator it,
                   std::vector<double> values) {
  double mean = 0.0;
  for (int count = 0; count < nr_values; ++count, ++it) {
    mean += values[*it];
  }
  mean /= nr_values;
  return mean;
}

double get_minimum(int nr_values, std::vector<int>::iterator it,
                   std::vector<double> values) {
  double min = values[*it++];
  ;
  for (int count = 1; count < nr_values; ++count, ++it) {
    double temp = values[*it];
    if (temp < min) {
      min = temp;
    }
  }
  return min;
}

void *_nearestNeighboursThreadEukledian(void *data) {
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
  struct timespec req = {0 /*req.tv_sec*/, 1L /*req.tv_nsec*/};
  GPSubData<int, double, double, double, int> *dat;
  dat = static_cast<GPSubData<int, double, double, double, int> *>(data);
  double *grdpnts = dat->GetData<0>();     // gridpoints that were split over threads
  double *all_grdpnts = dat->GetData<1>(); // the entire grid
  double *cutoff_in = dat->GetData<2>();
  double cutoff_used = (*cutoff_in) * (*cutoff_in); // eukledian dist only needs square

  int *config = dat->GetData<3>();
  int max_nr_neigh = *(config + 0);
  int nr_neigh = *(config + 1);
  bool sort_it = (*(config + 2) == 1);

  int *neigh_list = dat->GetDataOutput();
  const int nr_pnts = dat->GetNrOutput();
  const int nr_grdpnts = dat->GetNr<1>() / 3; // three doubles belong to one point
  const int progress = 250;
  const bool progress_reports = dat->GetProgressReports();
  int *progress_bar = dat->GetProgressBar();
  pthread_mutex_t *mut = dat->GetMutex();

  typedef std::pair<double, int> here_pair;
  here_pair *tbs; // tbs means "to be sorted" and will be sorted by distance
  tbs = (here_pair *)malloc(max_nr_neigh * sizeof(here_pair));

  double *p = grdpnts;
  int *nr = neigh_list;
  int *n = nr + 1;
  for (int ip = 0; ip < nr_pnts;) {
    for (int prog = 0; prog < progress && ip < nr_pnts; ++prog, ++ip) {
      double xpc = *(p + 0);
      double ypc = *(p + 1);
      double zpc = *(p + 2);
      double *c = all_grdpnts;
      int cnt = 0;
      std::pair<double, int> *t = tbs;
      for (int ic = 0; ic < nr_grdpnts && cnt < max_nr_neigh; ++ic, c += 3) {
        double dx = (*(c + 0) - xpc);
        double dy = (*(c + 1) - ypc);
        double dz = (*(c + 2) - zpc);
        double normsquared = dx * dx + dy * dy + dz * dz;
        if (normsquared < cutoff_used && normsquared > 0.0) {
          *t = std::make_pair(normsquared, ic);
          ++cnt;
          ++t;
        }
      }
      if (sort_it) {
        std::sort(tbs, tbs + cnt, sort_by_first);
      }
      *nr = cnt;
      t = tbs;
      for (int j = 0; j < cnt && j < nr_neigh; ++j, ++n, ++t) {
        *n = t->second;
      }
      for (int j = cnt; j < nr_neigh; ++j, ++n) {
        *n = -1;
      }
      nr = n;
      ++n;
      p += 3;
    }
    // the nanosleep function serves as a possible point where Ctrl-C
    // can interrupt the programme
    nanosleep(&req, (struct timespec *)NULL);
    if (progress_reports) {
      // the programme should not be cancelled while the mutex is locked
      pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);
      pthread_mutex_lock(mut);
      *progress_bar += progress;
      pthread_mutex_unlock(mut);
      pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    }
  }
  free(tbs);
  pthread_exit(NULL);
}

void *_nearestNeighboursThreadManhattanMultiple(void *data) {
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
  struct timespec req = {0 /*req.tv_sec*/, 1L /*req.tv_nsec*/};
  GPSubData<int, double, double, double, int> *dat;
  dat = static_cast<GPSubData<int, double, double, double, int> *>(data);
  double *grdpnts = dat->GetData<0>();     // gridpoints that were split over threads
  double *all_grdpnts = dat->GetData<1>(); // the entire grid

  double *cutoff_in = dat->GetData<2>();
  double cutoff_used[3];
  // manhattan distance uses three different values
  cutoff_used[0] = *(cutoff_in + 0);
  cutoff_used[1] = *(cutoff_in + 1);
  cutoff_used[2] = *(cutoff_in + 2);
  // this gets the absolute value faster than using abs
  cutoff_used[0] = cutoff_used[0] > 0 ? cutoff_used[0] : -cutoff_used[0];
  cutoff_used[1] = cutoff_used[1] > 0 ? cutoff_used[1] : -cutoff_used[1];
  cutoff_used[2] = cutoff_used[2] > 0 ? cutoff_used[2] : -cutoff_used[2];

  int *config = dat->GetData<3>();
  int max_nr_neigh = *(config + 0);
  int nr_neigh = *(config + 1);
  bool sort_it = (*(config + 2) == 1);

  int *neigh_list = dat->GetDataOutput();
  const int nr_pnts = dat->GetNrOutput();
  const int nr_grdpnts = dat->GetNr<1>() / 3; // three doubles belong to one point
  const int progress = 250;
  const bool progress_reports = dat->GetProgressReports();
  int *progress_bar = dat->GetProgressBar();
  pthread_mutex_t *mut = dat->GetMutex();

  typedef std::pair<double, int> here_pair;
  here_pair *tbs; // tbs means "to be sorted" and will be sorted by distance
  tbs = (here_pair *)malloc(max_nr_neigh * sizeof(here_pair));

  double *p = grdpnts;
  int *nr = neigh_list;
  int *n = nr + 1;
  for (int ip = 0; ip < nr_pnts;) {
    for (int prog = 0; prog < progress && ip < nr_pnts; ++prog, ++ip) {
      double xpc = *(p + 0);
      double ypc = *(p + 1);
      double zpc = *(p + 2);
      double *c = all_grdpnts;
      int cnt = 0;
      std::pair<double, int> *t = tbs;
      for (int ic = 0; ic < nr_grdpnts && cnt < max_nr_neigh; ++ic, c += 3) {
        double dx = (*(c + 0) - xpc);
        double dy = (*(c + 1) - ypc);
        double dz = (*(c + 2) - zpc);
        dx = dx > 0 ? dx : -dx;
        dy = dy > 0 ? dy : -dy;
        dz = dz > 0 ? dz : -dz;
        double norm = dx + dy + dz;
        if (dx < cutoff_used[0] && dy < cutoff_used[1] && dz < cutoff_used[2] &&
            norm > 0.0) {
          *t = std::make_pair(norm, ic);
          ++cnt;
          ++t;
        }
      }
      if (sort_it) {
        std::sort(tbs, tbs + cnt, sort_by_first);
      }
      *nr = cnt;
      t = tbs;
      for (int j = 0; j < cnt && j < nr_neigh; ++j, ++n, ++t) {
        *n = t->second;
      }
      for (int j = cnt; j < nr_neigh; ++j, ++n) {
        *n = -1;
      }
      nr = n;
      ++n;
      p += 3;
    }
    // the nanosleep function serves as a possible point where Ctrl-C
    // can interrupt the programme
    nanosleep(&req, (struct timespec *)NULL);
    if (progress_reports) {
      // the programme should not be cancelled while the mutex is locked
      pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);
      pthread_mutex_lock(mut);
      *progress_bar += progress;
      pthread_mutex_unlock(mut);
      pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    }
  }
  free(tbs);
  pthread_exit(NULL);
}

void *_nearestNeighboursThreadManhattanCombined(void *data) {
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
  struct timespec req = {0 /*req.tv_sec*/, 1L /*req.tv_nsec*/};
  GPSubData<int, double, double, double, int> *dat;
  dat = static_cast<GPSubData<int, double, double, double, int> *>(data);
  double *grdpnts = dat->GetData<0>();     // gridpoints that were split over threads
  double *all_grdpnts = dat->GetData<1>(); // the entire grid
  int *config = dat->GetData<3>();
  int max_nr_neigh = *(config + 0);
  int nr_neigh = *(config + 1);
  bool sort_it = (*(config + 2) == 1);
  double *cutoff_in = dat->GetData<2>();
  double cutoff_used = (*cutoff_in);
  cutoff_used = cutoff_used > 0 ? cutoff_used : -cutoff_used;
  int *neigh_list = dat->GetDataOutput();
  const int nr_pnts = dat->GetNrOutput();
  const int nr_grdpnts = dat->GetNr<1>() / 3; // three doubles belong to one point
  const int progress = 250;
  const bool progress_reports = dat->GetProgressReports();
  int *progress_bar = dat->GetProgressBar();
  pthread_mutex_t *mut = dat->GetMutex();

  typedef std::pair<double, int> here_pair;
  here_pair *tbs; // tbs means "to be sorted" and will be sorted by distance
  tbs = (here_pair *)malloc(max_nr_neigh * sizeof(here_pair));

  double *p = grdpnts;
  int *nr = neigh_list;
  int *n = nr + 1;
  for (int ip = 0; ip < nr_pnts;) {
    for (int prog = 0; prog < progress && ip < nr_pnts; ++prog, ++ip) {
      double xpc = *(p + 0);
      double ypc = *(p + 1);
      double zpc = *(p + 2);
      double *c = all_grdpnts;
      int cnt = 0;
      std::pair<double, int> *t = tbs;
      for (int ic = 0; ic < nr_grdpnts && cnt < max_nr_neigh; ++ic, c += 3) {
        double dx = (*(c + 0) - xpc);
        double dy = (*(c + 1) - ypc);
        double dz = (*(c + 2) - zpc);
        dx = dx > 0 ? dx : -dx;
        dy = dy > 0 ? dy : -dy;
        dz = dz > 0 ? dz : -dz;
        double norm = dx + dy + dz;
        if (norm < cutoff_used && norm > 0.0) {
          *t = std::make_pair(norm, ic);
          ++cnt;
          ++t;
        }
      }
      if (sort_it) {
        std::sort(tbs, tbs + cnt, sort_by_first);
      }
      *nr = cnt;
      t = tbs;
      for (int j = 0; j < cnt && j < nr_neigh; ++j, ++n, ++t) {
        *n = t->second;
      }
      for (int j = cnt; j < nr_neigh; ++j, ++n) {
        *n = -1;
      }
      nr = n;
      ++n;
      p += 3;
    }
    // the nanosleep function serves as a possible point where Ctrl-C
    // can interrupt the programme
    nanosleep(&req, (struct timespec *)NULL);
    if (progress_reports) {
      // the programme should not be cancelled while the mutex is locked
      pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);
      pthread_mutex_lock(mut);
      *progress_bar += progress;
      pthread_mutex_unlock(mut);
      pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    }
  }
  free(tbs);
  pthread_exit(NULL);
}

void make_neighbour_list_irregular(bool progress_reports, int nr_gridpoints,
                                   int max_nr_neighbours, int nr_neighbours,
                                   int cutoff_type, std::vector<double> points,
                                   std::vector<double> distance_cutoff,
                                   std::vector<int> *neighbour_list, bool sort_it) {
  // initialize everything
  PG globals;
  init_parallel_generic(&progress_reports, &globals);

  std::vector<double> vec_cutoff;
  void *(*neighbour_func)(void *);

  switch (cutoff_type) {
  case 1:
    neighbour_func = _nearestNeighboursThreadEukledian;
    vec_cutoff.reserve(1);
    vec_cutoff.push_back(distance_cutoff.at(0));
    break;
  case 2:
    neighbour_func = _nearestNeighboursThreadManhattanMultiple;
    vec_cutoff.reserve(3);
    vec_cutoff.push_back(distance_cutoff.at(0));
    vec_cutoff.push_back(distance_cutoff.at(1));
    vec_cutoff.push_back(distance_cutoff.at(2));
    break;
  case 3:
    neighbour_func = _nearestNeighboursThreadManhattanCombined;
    vec_cutoff.reserve(1);
    vec_cutoff.push_back(distance_cutoff.at(0));
    break;
  default:
    throw std::invalid_argument(
        "Supported metrics: 1: Eukledian, 2: Manhattan in multiple directions, 3: "
        "Manhattan combined. Choose one of them.");
    break;
  }

  std::vector<int> vec_config;
  vec_config.reserve(3);
  vec_config.push_back(max_nr_neighbours);
  vec_config.push_back(nr_neighbours);
  vec_config.push_back(sort_it ? 1 : 0);

  // reserve data structures and fill them with input
  tuple_of_vectors<double, double, double, int> input;
  input = std::make_tuple(points, points, vec_cutoff, vec_config);

  const int split_factor_in = 3;
  const int split_factor_out = nr_neighbours + 1;

  // fill class that holds data for each thread
  GPData<int, double, double, double, int> *data;
  try {
    data = new GPData<int, double, double, double, int>(progress_reports,
                                                        globals.nr_threads,
                                                        input,
                                                        neighbour_list,
                                                        &(globals.mutex),
                                                        &(globals.progress_bar),
                                                        split_factor_in,
                                                        split_factor_out,
                                                        true);
  } catch (const std::invalid_argument &e) {
    throw;
  }
  // perform computation
  do_parallel_generic<int, double, double, double, int>(
      neighbour_func, &globals, progress_reports, nr_gridpoints, data);
  // transfer output data
  data->TransferOutput();
  // clean up
  delete data;
  finalize_parallel_generic(progress_reports, &globals);
  pg_global = NULL;
}

/**
 * \brief A class that aides in finding indices of neighbours to a point on a regular
 * grid.
 *
 * The class is initialized using the grids dimensions (nx, ny, nz: number of
 * points in each direction).  Then, it is passed the one-dimensional index
 * of a point (starting at 0 and ending at nx*ny*nz-1). Then, when given a
 * displacement (in the form of index offsets in the 3 Cartesian directions) it
 * returns the one-dimensional index of that point (if it exists in the grid).
 * I implemented it this way because one-dimensional indices have to be used
 * with flat data structures (which are easier to handle, IMHO) but it is easier
 * to think in terms of three-dimensional indices when it comes to regular grids.
 */
class Slices {
private:
  /*! int - number of points in x-direction */
  const int _nx;
  /*! int - number of points in y-direction */
  const int _ny;
  /*! int - number of points in z-direction */
  const int _nz;
  /*! int - total number of points in the grid */
  const int _nr;
  /*! int - x-part of the 3d index of a given point */
  int _sx;
  /*! int - y-part of the 3d index of a given point */
  int _sy;
  /*! int - z-part of the 3d index of a given point */
  int _sz;

public:
  /** \brief Constructor
   * \param nx int - number of points in x-direction
   * \param ny int - number of points in y-direction
   * \param nz int - number of points in z-direction
   */
  Slices(int nx, int ny, int nz)
      : _nx(nx), _ny(ny), _nz(nz), _nr(nx * ny * nz), _sx(-1), _sy(-1), _sz(-1) {}

  /**
   * \brief Declare a reference point.
   *
   * When passing a 3d displacement to \a GetNeighbourIndex, the displacements
   * are taken relative to the point declared in this function.
   *
   * \param index int - one dimensional index of the point
   * \return whether or not the poin is on the grid
   */
  bool SetPoint(int index) {
    if (index < 0 || index > _nr) {
      return false;
    }
    _sz = index % _nz;
    _sy = (index / _nz) % _ny;
    _sx = ((index / _nz) / _ny) % _nx;
    if (_sx < _nx && _sy < _ny && _sz < _nz) {
      return true;
    } else {
      return false;
    }
  }

  /**
   * \brief Get the one dimensional index of a point relative to a central point.
   *
   * The central point is declared using \a SetPoint.
   *
   * \param dx int - index displacement in x-direction
   * \param dy int - index displacement in y-direction
   * \param dz int - index displacement in z-direction
   * \return the one dimensional index of the point
   */
  int GetNeighbourIndex(int dx, int dy, int dz) {
    // a point is never its own neighbour
    if (dx == 0 && dy == 0 && dz == 0) {
      return -1;
    }
    int neigh_x, neigh_y, neigh_z;
    neigh_x = _sx + dx;
    if (neigh_x < 0 || neigh_x >= _nx) {
      return -1;
    }
    neigh_y = _sy + dy;
    if (neigh_y < 0 || neigh_y >= _ny) {
      return -1;
    }
    neigh_z = _sz + dz;
    if (neigh_z < 0 || neigh_z >= _nz) {
      return -1;
    }
    // return the actual index of the neighbour
    return (neigh_x * _ny + neigh_y) * _nz + neigh_z;
  }
};

void make_neighbour_list_regular(bool progress_reports, bool exclude_border,
                                 int nr_gridpoints_x, int nr_gridpoints_y,
                                 int nr_gridpoints_z, int nr_neighbour_shells,
                                 std::vector<int> *neighbour_list) {
  // x is slow, y is middle and z is fast index
  // only one thread used because this is fast

  // initialize everything
  if (nr_gridpoints_x < 0) {
    fprintf(stderr, "ERROR: nr_gridpoints_x must be non-negative.");
    return;
  }
  if (nr_gridpoints_y < 0) {
    fprintf(stderr, "ERROR: nr_gridpoints_y must be non-negative.");
    return;
  }
  if (nr_gridpoints_z < 0) {
    fprintf(stderr, "ERROR: nr_gridpoints_z must be non-negative.");
    return;
  }
  if (nr_neighbour_shells < 0) {
    fprintf(stderr, "ERROR: nr_neighbour_shells must be non-negative.");
    return;
  }

  bool too_many = false;
  if (nr_gridpoints_x > std::numeric_limits<int>::max() / nr_gridpoints_y) {
    too_many = true;
  }
  if (too_many || (nr_gridpoints_x * nr_gridpoints_y >
                   std::numeric_limits<int>::max() / nr_gridpoints_z)) {
    too_many = true;
  }
  if (too_many) {
    fprintf(stderr, "ERROR: too many points in grid, cannot index them with integers.");
    return;
  }

  int ng =
      nr_gridpoints_x * nr_gridpoints_y * nr_gridpoints_z; // number of total gridpoints

  // determine how many neighbours are there at max in all shells. Compute as:
  // nr_neigh_in_shells = ((2*(nr_neighbour_shells)+1)^3)-1, i.e., as a cube with the
  // central point removed
  unsigned int ni = (2 * nr_neighbour_shells) + 1; // nr_neigh_in_shells
  ni *= ni * ni;                                   // ni^3
  --ni; // before this step, ni will always be >=1

  neighbour_list->reserve(ng * (ni + 1));
  neighbour_list->assign(ng * (ni + 1), -1); // assign -1 to all elements

  std::vector<int>::iterator startit;
  startit = neighbour_list->begin();

  // this data structure is used to get the x, y and z indices from the overall index
  Slices slices = Slices(nr_gridpoints_x, nr_gridpoints_y, nr_gridpoints_z);

  if (progress_reports) {
    printf("Starting generation of neighbour list of regular grid.\n");
    printf("Prog: %.2f\n", 0.0);
    fflush(stdout);
  }
  // loop over all points in the grid
  for (int p = 0; p < ng; ++p, startit += (ni + 1)) {
    // save the data for the current point to this vector
    std::vector<int> point_data;
    point_data.reserve(ni + 1);
    point_data.assign(ni + 1, -1);
    // start at the first element for this point
    std::vector<int>::iterator point_it = point_data.begin() + 1;
    // find the number of the slice in each direction, i.e., the x, y and z index from
    // the overall one
    if (!slices.SetPoint(p)) {
      fprintf(stderr, "ERROR treating point number %d, skipping.", p);
      continue;
    }
    // this will count how many neighbours already have been found
    unsigned int nr_neigh = 0;
    // check whether has neighbour in all directions
    for (int ix = -nr_neighbour_shells; ix <= nr_neighbour_shells; ++ix) {
      for (int iy = -nr_neighbour_shells; iy <= nr_neighbour_shells; ++iy) {
        for (int iz = -nr_neighbour_shells; iz <= nr_neighbour_shells; ++iz) {
          int index;
          if ((index = slices.GetNeighbourIndex(ix, iy, iz)) != -1) {
            ++nr_neigh;
            *point_it = index;
            ++point_it;
          }
        } // iz
      }   // iy
    }     // ix
    point_it = point_data.begin();
    *point_it = (int)nr_neigh;
    if (!exclude_border || nr_neigh == ni) {
      // assign the actual number of neighbours to the appropriate place and copy data
      for (std::vector<int>::iterator it = startit; point_it != point_data.end();
           ++it, ++point_it) {
        *it = *point_it;
      }
    } else {
      *startit = 0;
    }
    if (progress_reports) {
      if (p % 500 == 0) {
        fprintf(stdout, "%c[2K\r", 27);
        printf("Prog: %.2f", 100.0 * p / ng);
        fflush(stdout);
      }
    }
  }
  if (progress_reports) {
    fprintf(stdout, "%c[2K\r", 27);
    printf("Prog: %.2f\n", 100.0);
    fflush(stdout);
  }
}

void local_minima_from_neighbour_list(bool progress_reports, int nr_neighbours,
                                      int nr_values, std::vector<int> neighbour_list,
                                      std::vector<double> values,
                                      std::vector<int> *minima,
                                      std::vector<double> degeneration_cutoffs,
                                      bool use_upper_cutoff, bool use_lower_cutoff,
                                      double upper_cutoff, double lower_cutoff,
                                      int sort_it, std::vector<double> *depths) {

  // the iterator over the neighbour list
  // will be kept until minima have been determined
  std::vector<int>::iterator neigh_it;
  std::vector<int>::iterator next_neigh_it = neighbour_list.begin();

  // the iterator over the volumetric data values
  std::vector<double>::iterator val_it;

  // this will contain "true" at the corresponding index
  // if a point is a local minimum
  std::vector<bool> is_min_vec;
  is_min_vec.reserve(nr_values);

  int val_count = 0;

  if (progress_reports) {
    printf("Prog: %.2f", 0.0);
    fflush(stdout);
  }
  for (val_it = values.begin(); val_it != values.end(); ++val_it, ++val_count) {
    bool is_min =
        true; // innocent until proven guilty. the americans could learn from that ^^.
    int nr_neigh = *next_neigh_it;
    std::vector<int>::iterator temp_it = next_neigh_it + 1;
    neigh_it = next_neigh_it + 1;
    next_neigh_it += nr_neighbours + 1;
    double here_value = *val_it;
    bool screened;
    screened = use_upper_cutoff ? here_value > upper_cutoff : false;
    screened = use_lower_cutoff ? (screened || here_value < lower_cutoff) : screened;
    if (screened || nr_neigh <= 0) {
      is_min = false;
    } else {
      here_value += degeneration_cutoffs[0];
      for (; nr_neigh > 0; --nr_neigh, ++neigh_it) {
        double neigh_value = values[*neigh_it];
        bool neigh_screened;
        // I do not want to consider points whose neighbours are screened to
        // be cancidates for minima because they might not be true but artificial
        // minima
        neigh_screened = use_upper_cutoff ? neigh_value > upper_cutoff : false;
        neigh_screened = use_lower_cutoff
                             ? (neigh_screened || neigh_value < lower_cutoff)
                             : neigh_screened;
        if (neigh_screened || (here_value > neigh_value)) {
          is_min = false;
          break;
        }
      }
    }
    is_min_vec.push_back(is_min);
    if (progress_reports) {
      if (val_count % 500 == 0) {
        fprintf(stdout, "%c[2K\r", 27);
        printf("Prog: %.2f", 100.0 * val_count / nr_values);
        fflush(stdout);
      }
    }
  }
  if (progress_reports) {
    fprintf(stdout, "%c[2K\r", 27);
    printf("Prog: %.2f\n", 100.0);
    fflush(stdout);
  }
  int nr_minima = std::count_if(
      is_min_vec.begin(), is_min_vec.end(), [](bool b) -> bool { return b; });

  minima->reserve(nr_minima);

  int is_min_count = 0;
  for (std::vector<bool>::iterator bool_it = is_min_vec.begin();
       bool_it != is_min_vec.end();
       ++bool_it, ++is_min_count) {
    if (*bool_it) {
      minima->push_back(is_min_count);
    }
  }

  if (sort_it != 0) {
    double (*get_compare_value_func)(
        int, std::vector<int>::iterator, std::vector<double>);
    switch (sort_it) {
    case 1:
      get_compare_value_func = get_average;
      break;
    case 2:
      get_compare_value_func = get_minimum;
      break;
    default:
      throw std::invalid_argument("Supported comparison algorithms: 1: average, 2: "
                                  "lowest value of all neighbours. Choose one.");
    }
    std::vector<std::pair<double, int>> tbs;
    tbs.reserve(nr_minima);
    for (std::vector<int>::iterator int_it = minima->begin(); int_it != minima->end();
         ++int_it) {
      double compare_value, here_value;
      // this points to the start of data associated with the (*int_it)-th point
      neigh_it = neighbour_list.begin() + (*int_it * (nr_neighbours + 1));
      compare_value = get_compare_value_func(*neigh_it, neigh_it + 1, values);
      here_value = values[*int_it];
      tbs.push_back(std::make_pair(compare_value - here_value, *int_it));
    }
    std::sort(tbs.begin(), tbs.end(), sort_by_first);
    minima->clear();
    minima->reserve(nr_minima);
    // since the sorting happens in ascending order but the result shall have descending
    // order, iterate through the vector in reverse order
    for (std::vector<std::pair<double, int>>::reverse_iterator pair_it = tbs.rbegin();
         pair_it != tbs.rend();
         ++pair_it) {
      minima->push_back(pair_it->second);
    }

    // if requested, return the depths that were obtained
    if (depths) {
      depths->reserve(nr_minima);
      for (std::vector<std::pair<double, int>>::reverse_iterator pair_it = tbs.rbegin();
           pair_it != tbs.rend();
           ++pair_it) {
        depths->push_back(pair_it->first);
      }
    }
  }
}
