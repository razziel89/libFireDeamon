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
#include <FireDeamon/core/constants.h>
#include <FireDeamon/core/electron_density.h>
#include <FireDeamon/core/parallel_generic.h>
#include <math.h>
#include <pthread.h>
#include <stdexcept>
#include <time.h>
#include <tuple>
#include <vector>

// const static double pi_to_3_half = 5.5683279968317078453;
// const static double sqrt_pihalf_to_3_4 = 1.4031041455342160267;
// const static double two_div_by_pi_to_three_fourth = 0.71270547035499016035;
// const static double sqrt2 = 1.41421356237309504880;

inline double power(double base, unsigned int exp) {
  double result = 1.0;
  while (exp) {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
}

void *_electronDensityThreadCutoff(void *data) {
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
  struct timespec req = {0 /*req.tv_sec*/, 1L /*req.tv_nsec*/};
  // req.tv_sec = 0;
  // req.tv_nsec = 1L;
  GPSubData<double, double, double, int, double, double, double, double> *dat;
  dat = static_cast<
      GPSubData<double, double, double, int, double, double, double, double> *>(data);

  double *grdpnts = dat->GetData<0>();  // gridpoints
  double *prm_cent = dat->GetData<1>(); // centers of primitives
  int *prm_ang = dat->GetData<2>();     // angular exponents of primitives
  double *prm_exp = dat->GetData<3>();  // exponents of primitives
  double *prm_coef = dat->GetData<4>(); // contraction coefficients of primitives
  double *mo_coef = dat->GetData<5>();  // molecular orbital coefficients
  double *cutoff = dat->GetData<6>(); // cutoff above which no density will be computed
  double cutoff_2 = (*cutoff) * (*cutoff); // only the square is needed
  double *dens = dat->GetDataOutput();

  const int nr_coef = dat->GetNr<4>();
  const int nr_dens = dat->GetNrOutput();
  const int progress = 5000;
  const bool progress_reports = dat->GetProgressReports();
  int *progress_bar = dat->GetProgressBar();
  pthread_mutex_t *mut = dat->GetMutex();

  double *gp = grdpnts;
  double *density = dens;
  for (int dp = 0; dp < nr_dens;) {
    for (int prog = 0; prog < progress && dp < nr_dens; ++prog, ++dp) {
      double sum = 0.0;
      double xgp = *gp;
      ++gp;
      double ygp = *gp;
      ++gp;
      double zgp = *gp;
      ++gp;
      double *pcent = prm_cent;
      int *pang = prm_ang;
      double *pexp = prm_exp;
      double *pcoef = prm_coef;
      double *mo = mo_coef;
      for (int ip = 0; ip < nr_coef; ++ip) {
        double dx = xgp - (*(pcent + 0));
        double dy = ygp - (*(pcent + 1));
        double dz = zgp - (*(pcent + 2));
        double dist_2 = (dx * dx + dy * dy + dz * dz);
        if (dist_2 < cutoff_2) {
          double exp_factor = exp(-(*pexp) * dist_2);
          double prefac;
          prefac = power(dx, *(pang + 0));
          prefac *= power(dy, *(pang + 1));
          prefac *= power(dz, *(pang + 2));
          sum += (*pcoef) * (*mo) * prefac * exp_factor;
        }
        ++pcoef;
        ++mo;
        ++pexp;
        pang += 3;
        pcent += 3;
      }
      *density = sum * sum;
      ++density;
    }
    // the nanosleep function serves as a possible point where Ctrl-C
    // can interrupt the programme
    nanosleep(&req, (struct timespec *)NULL);
    // the programme should not be cancelled while the mutex is locked
    if (progress_reports) {
      pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);
      pthread_mutex_lock(mut);
      *progress_bar += progress;
      pthread_mutex_unlock(mut);
      pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    }
  }
  pthread_exit(NULL);
}

void *_electronDensityThreadNoCutoff(void *data) {
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
  struct timespec req = {0 /*req.tv_sec*/, 1L /*req.tv_nsec*/};
  // req.tv_sec = 0;
  // req.tv_nsec = 1L;
  GPSubData<double, double, double, int, double, double, double, double> *dat;
  dat = static_cast<
      GPSubData<double, double, double, int, double, double, double, double> *>(data);

  double *grdpnts = dat->GetData<0>();  // gridpoints
  double *prm_cent = dat->GetData<1>(); // centers of primitives
  int *prm_ang = dat->GetData<2>();     // angular exponents of primitives
  double *prm_exp = dat->GetData<3>();  // exponents of primitives
  double *prm_coef = dat->GetData<4>(); // contraction coefficients of primitives
  double *mo_coef = dat->GetData<5>();  // molecular orbital coefficients
  double *dens = dat->GetDataOutput();

  const int nr_coef = dat->GetNr<4>();
  const int nr_dens = dat->GetNrOutput();
  // const int sub_nr = dat->GetSubNr();
  const int progress = 5000;
  const bool progress_reports = dat->GetProgressReports();
  int *progress_bar = dat->GetProgressBar();
  pthread_mutex_t *mut = dat->GetMutex();

  double *gp = grdpnts;
  double *density = dens;
  for (int dp = 0; dp < nr_dens;) {
    for (int prog = 0; prog < progress && dp < nr_dens; ++prog, ++dp) {
      double sum = 0.0;
      double xgp = *gp;
      ++gp;
      double ygp = *gp;
      ++gp;
      double zgp = *gp;
      ++gp;
      double *pcent = prm_cent;
      int *pang = prm_ang;
      double *pexp = prm_exp;
      double *pcoef = prm_coef;
      double *mo = mo_coef;
      for (int ip = 0; ip < nr_coef; ++ip) {
        double dx = xgp - (*pcent);
        ++pcent;
        double dy = ygp - (*pcent);
        ++pcent;
        double dz = zgp - (*pcent);
        ++pcent;
        double dist_2 = (dx * dx + dy * dy + dz * dz);
        double exp_factor = exp(-(*pexp) * dist_2);
        double prefac;
        prefac = power(dx, *(pang + 0));
        prefac *= power(dy, *(pang + 1));
        prefac *= power(dz, *(pang + 2));
        sum += (*pcoef) * (*mo) * prefac * exp_factor;
        ++pcoef;
        ++mo;
        ++pexp;
        pang += 3;
      }
      *density = sum * sum;
      ++density;
    }
    // the nanosleep function serves as a possible point where Ctrl-C
    // can interrupt the programme
    nanosleep(&req, (struct timespec *)NULL);
    // the programme should not be cancelled while the mutex is locked
    if (progress_reports) {
      pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);
      pthread_mutex_lock(mut);
      *progress_bar += progress;
      pthread_mutex_unlock(mut);
      pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    }
  }
  pthread_exit(NULL);
}

////this name is a short form of: one_div_by_sqrt_double_factorial_of_2aminus1
//// which is 1.0/((2*a-1)!!)
// double odbsdfo2[] = {
//    1.0/sqrt(1.0),
//    1.0/sqrt(1.0),
//    1.0/sqrt(3.0),
//    1.0/sqrt(15.0),
//    1.0/sqrt(105.0),
//    1.0/sqrt(945.0),
//    1.0/sqrt(10395.0),
//    1.0/sqrt(135135.0),
//    1.0/sqrt(2027025.0),
//    1.0/sqrt(34459425.0),
//    1.0/sqrt(654729075.0),
//    1.0/sqrt(13749310575.0),
//    1.0/sqrt(316234143225.0),
//    1.0/sqrt(7905853580625.0),
//    1.0/sqrt(213458046676875.0),
//    1.0/sqrt(6190283353629375.0),
//    1.0/sqrt(191898783962510625.0),
//    1.0/sqrt(6332659870762850625.0)
//};

void normalize_gaussians(std::vector<double> *prefactor, std::vector<double> exponent,
                         std::vector<int> angular) {
  if (prefactor->size() != exponent.size()) {
    throw std::invalid_argument("Lengths of prefactor, exponent vectors do not equal.");
  }
  if (prefactor->size() * 3 != angular.size()) {
    throw std::invalid_argument(
        "Lengths of prefactor (times 3) and angular vectors do not equal.");
  }
  std::vector<double>::iterator preit, expit;
  std::vector<int>::iterator angit;
  for (preit = prefactor->begin(), expit = exponent.begin(), angit = angular.begin();
       preit != prefactor->end();
       ++preit, ++expit, angit += 3) {
    // printf("%f ",*preit); fflush(stdout);
    //*preit *= two_div_by_pi_to_three_fourth * pow(*expit,0.75);
    //*preit *= two_div_by_pi_to_three_fourth * sqrt2 * pow(*expit,1.5);
    // printf("%f\n",*preit); fflush(stdout);
    int ax = *(angit + 0);
    int ay = *(angit + 1);
    int az = *(angit + 2);
    // odbsdfo2 is short for one_div_by_sqrt_double_factorial_of_2aminus1 (see above)
    *preit *= two_div_by_pi_to_three_fourth * pow(*expit, 0.75) *
              pow(2.0 * 2.0 * (*expit), 0.5 * (ax + ay + az)) * odbsdfo2[ax] *
              odbsdfo2[ay] * odbsdfo2[az];
  }
}

void electron_density(bool progress_reports, int num_gridpoints,
                      std::vector<double> prim_centers,
                      std::vector<double> prim_exponents,
                      std::vector<double> prim_coefficients,
                      std::vector<int> prim_angular, std::vector<double> density_grid,
                      std::vector<double> mo_coefficients, std::vector<double> *density,
                      double cutoff) {
  // initialize everything
  PG globals;
  init_parallel_generic(&progress_reports, &globals);

  // check whether a cutoff shall be considered
  bool use_cutoff = (cutoff > 0.0);
  std::vector<double> vec_cutoff;
  vec_cutoff.push_back(cutoff);

  const int split_factor = 3;

  // reserve data structures and fill them with input
  tuple_of_vectors<double, double, int, double, double, double, double> input;
  input = std::make_tuple(density_grid,
                          prim_centers,
                          prim_angular,
                          prim_exponents,
                          prim_coefficients,
                          mo_coefficients,
                          vec_cutoff);

  // fill class that holds data for each thread
  GPData<double, double, double, int, double, double, double, double> *data;
  try {
    data = new GPData<double, double, double, int, double, double, double, double>(
        progress_reports,
        globals.nr_threads,
        input,
        density,
        &(globals.mutex),
        &(globals.progress_bar),
        split_factor,
        1,
        /*interlace=*/true);
  } catch (const std::invalid_argument &e) {
    throw;
  }
  fflush(stdout);
  // perform computation
  if (use_cutoff) {
    do_parallel_generic<double, double, double, int, double, double, double, double>(
        _electronDensityThreadCutoff, &globals, progress_reports, num_gridpoints, data);
  } else {
    do_parallel_generic<double, double, double, int, double, double, double, double>(
        _electronDensityThreadNoCutoff,
        &globals,
        progress_reports,
        num_gridpoints,
        data);
  }
  // transfer output data
  data->TransferOutput();
  // clean up
  delete data;
  finalize_parallel_generic(progress_reports, &globals);
  pg_global = NULL;
}
