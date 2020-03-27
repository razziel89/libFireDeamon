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
#include <FireDeamon/core/irregular_grid_interpolation.h>
#include <FireDeamon/core/parallel_generic.h>
#include <limits>
#include <math.h>
#include <pthread.h>
#include <stdexcept>
#include <time.h>
#include <tuple>
#include <vector>

void *_nearestInterpolationThread(void *data) {
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
  struct timespec req = {0 /*req.tv_sec*/, 1L /*req.tv_nsec*/};
  GPSubData<double, double, double, double, double, int> *dat;
  dat = static_cast<GPSubData<double, double, double, double, double, int> *>(data);
  double *inpnts = dat->GetData<0>(); // interpolation points
  double *pnts = dat->GetData<1>();
  double *vals = dat->GetData<2>();
  double *interp = dat->GetDataOutput();
  const int nr_pnts = dat->GetNr<1>();
  const int nr_interp = dat->GetNrOutput();
  const int progress = 250;
  const bool progress_reports = dat->GetProgressReports();
  int *progress_bar = dat->GetProgressBar();
  pthread_mutex_t *mut = dat->GetMutex();

  double *p = interp;
  double *pc = inpnts;
  for (int ip = 0; ip < nr_interp;) {
    for (int prog = 0; prog < progress && ip < nr_interp; ++prog, ++ip, ++p) {
      double min_dist_squ = std::numeric_limits<double>::max();
      double min_dist_val = 0.0;
      double xpc = *(pc + 0);
      double ypc = *(pc + 1);
      double zpc = *(pc + 2);
      double *c = pnts;
      double *v = vals;
      for (int ic = 0; ic < nr_pnts; ic += 3) {
        double dx = (*(c + 0) - xpc);
        double dy = (*(c + 1) - ypc);
        double dz = (*(c + 2) - zpc);
        double normsquared = dx * dx + dy * dy + dz * dz;
        if (normsquared < min_dist_squ) {
          min_dist_squ = normsquared;
          min_dist_val = *v;
        }
        ++v;
        c += 3;
      }
      *p = min_dist_val;
      pc += 3;
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
  pthread_exit(NULL);
}

void *_inverseDistanceWeightingInterpolationNoCutoffThread(void *data) {
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
  struct timespec req = {0 /*req.tv_sec*/, 1L /*req.tv_nsec*/};
  GPSubData<double, double, double, double, double, int> *dat;
  dat = static_cast<GPSubData<double, double, double, double, double, int> *>(data);
  double *inpnts = dat->GetData<0>(); // interpolation points
  double *pnts = dat->GetData<1>();
  double *vals = dat->GetData<2>();
  int *config = dat->GetData<4>();
  const int distance_exponent = config[0];
  const int distance_function = config[1];
  const double combined_exponent = -distance_exponent * 1.0 / distance_function;
  double *interp = dat->GetDataOutput();
  const int nr_pnts = dat->GetNr<1>();
  const int nr_interp = dat->GetNrOutput();
  const int progress = 25;
  const bool progress_reports = dat->GetProgressReports();
  int *progress_bar = dat->GetProgressBar();
  pthread_mutex_t *mut = dat->GetMutex();

  double *p = interp;
  double *pc = inpnts;
  for (int ip = 0; ip < nr_interp;) {
    for (int prog = 0; prog < progress && ip < nr_interp; ++prog, ++ip, ++p) {
      long double sum_weights = 0.0L;
      long double sum_weiths_times_value = 0.0L;
      double xpc = *(pc + 0);
      double ypc = *(pc + 1);
      double zpc = *(pc + 2);
      double *c = pnts;
      double *v = vals;
      for (int ic = 0; ic < nr_pnts; ic += 3) {
        double dx = *(c + 0) - xpc;
        double dy = *(c + 1) - ypc;
        double dz = *(c + 2) - zpc;
        double tx, ty, tz, ttx, tty, ttz;
        tx = ((dx > 0) ? dx : -dx);
        ty = ((dy > 0) ? dy : -dy);
        tz = ((dz > 0) ? dz : -dz);
        ttx = tx;
        tty = ty;
        ttz = tz;
        for (int j = 1; j < distance_function; ++j) {
          ttx *= tx;
          tty *= ty;
          ttz *= tz;
        }
        double weight = pow(ttx + tty + ttz, combined_exponent);
        sum_weights += (long double)weight;
        sum_weiths_times_value += (long double)(weight * (*v));
        ++v;
        c += 3;
      }
      *p = (double)(sum_weiths_times_value / sum_weights);
      pc += 3;
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
  pthread_exit(NULL);
}

void *_inverseDistanceWeightingInterpolationCutoffThread(void *data) {
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
  struct timespec req = {0 /*req.tv_sec*/, 1L /*req.tv_nsec*/};
  GPSubData<double, double, double, double, double, int> *dat;
  dat = static_cast<GPSubData<double, double, double, double, double, int> *>(data);
  double *inpnts = dat->GetData<0>(); // interpolation points
  double *pnts = dat->GetData<1>();
  double *vals = dat->GetData<2>();
  double *cutoff = dat->GetData<3>();
  int *config = dat->GetData<4>();
  double cutoff_2 = (*cutoff) * (*cutoff);
  const int distance_exponent = config[0];
  const int distance_function = config[1];
  const double combined_exponent = -distance_exponent * 1.0 / distance_function;
  double *interp = dat->GetDataOutput();
  const int nr_pnts = dat->GetNr<1>();
  const int nr_interp = dat->GetNrOutput();
  const int progress = 25;
  const bool progress_reports = dat->GetProgressReports();
  int *progress_bar = dat->GetProgressBar();
  pthread_mutex_t *mut = dat->GetMutex();

  double *p = interp;
  double *pc = inpnts;
  for (int ip = 0; ip < nr_interp;) {
    for (int prog = 0; prog < progress && ip < nr_interp; ++prog, ++ip, ++p) {
      long double sum_weights = 0.0L;
      long double sum_weiths_times_value = 0.0L;
      double xpc = *(pc + 0);
      double ypc = *(pc + 1);
      double zpc = *(pc + 2);
      double *c = pnts;
      double *v = vals;
      for (int ic = 0; ic < nr_pnts; ic += 3) {
        double dx = *(c + 0) - xpc;
        double dy = *(c + 1) - ypc;
        double dz = *(c + 2) - zpc;
        if (dx * dx + dy * dy + dz * dz < cutoff_2) {
          double tx, ty, tz, ttx, tty, ttz;
          tx = ((dx > 0) ? dx : -dx);
          ty = ((dy > 0) ? dy : -dy);
          tz = ((dz > 0) ? dz : -dz);
          ttx = tx;
          tty = ty;
          ttz = tz;
          for (int j = 1; j < distance_function; ++j) {
            ttx *= tx;
            tty *= ty;
            ttz *= tz;
          }
          double weight = pow(ttx + tty + ttz, combined_exponent);
          sum_weights += (long double)weight;
          sum_weiths_times_value += (long double)(weight * (*v));
        }
        ++v;
        c += 3;
      }
      *p = (double)(sum_weiths_times_value / sum_weights);
      pc += 3;
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
  pthread_exit(NULL);
}

void generic_interpolation(bool progress_reports, int num_interpolation_points,
                           std::vector<double> points, std::vector<double> values,
                           std::vector<double> interpolation_points,
                           std::vector<double> *interpolation, int interpolation_type,
                           int distance_exponent, int distance_function,
                           double cutoff) {
  void *(*interp_func)(void *);
  // initialize everything
  PG globals;
  init_parallel_generic(&progress_reports, &globals);

  // reserve data structures and fill them with input
  const int split_factor = 3;

  // check whether a cutoff shall be considered
  bool use_cutoff = (cutoff > 0.0);
  std::vector<double> vec_cutoff;
  vec_cutoff.push_back(cutoff);

  std::vector<int> config;
  if (interpolation_type == 2) {
    config.reserve(2);
    config.push_back(distance_exponent);
    config.push_back(distance_function);
  } else {
    config.reserve(0);
  }

  tuple_of_vectors<double, double, double, double, int> input;
  input = std::make_tuple(interpolation_points, points, values, vec_cutoff, config);

  // fill class that holds data for each thread
  GPData<double, double, double, double, double, int> *data;
  switch (interpolation_type) {
  case 1:
    if (progress_reports) {
      fprintf(stdout, "Using nearest-neighbour interpolation.\n");
    }
    interp_func = _nearestInterpolationThread;
    try {
      data = new GPData<double, double, double, double, double, int>(
          progress_reports,
          globals.nr_threads,
          input,
          interpolation,
          &(globals.mutex),
          &(globals.progress_bar),
          split_factor,
          1,
          false);
    } catch (const std::invalid_argument &e) {
      throw;
    }
    break;
  case 2:
    assert(distance_exponent >= 0.0 && distance_function >= 1);
    if (progress_reports) {
      fprintf(stdout,
              "Using inverse distance weighting with p=%d and the ",
              distance_exponent);
      switch (distance_function % 10) {
      case 1:
        fprintf(stdout, "%dst", distance_function);
        break;
      case 2:
        fprintf(stdout, "%dnd", distance_function);
        break;
      case 3:
        fprintf(stdout, "%drd", distance_function);
        break;
      default:
        fprintf(stdout, "%dth", distance_function);
        break;
      }
      fprintf(stdout, " root.\n");
    }
    if (use_cutoff) {
      interp_func = _inverseDistanceWeightingInterpolationCutoffThread;
    } else {
      interp_func = _inverseDistanceWeightingInterpolationNoCutoffThread;
    }
    try {
      data = new GPData<double, double, double, double, double, int>(
          progress_reports,
          globals.nr_threads,
          input,
          interpolation,
          &(globals.mutex),
          &(globals.progress_bar),
          split_factor,
          1,
          false);
    } catch (const std::invalid_argument &e) {
      throw;
    }
    break;
  default:
    throw std::invalid_argument("Supported interpolation methods: 1: nearest "
                                "neighbour, 2: inverse distance. Choose one of them.");
    break;
  }
  fflush(stdout);
  // perform computation
  do_parallel_generic<double, double, double, double, double, int>(
      interp_func, &globals, progress_reports, num_interpolation_points, data);
  // transfer output data
  data->TransferOutput();
  // clean up
  delete data;
  finalize_parallel_generic(progress_reports, &globals);
  pg_global = NULL;
}
