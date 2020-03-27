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
#include <FireDeamon/core/electrostatic_potential_orbitals.h>
#include <FireDeamon/core/halfnum/angular_integral.h>
#include <FireDeamon/core/halfnum/radial_integral.h>
#include <FireDeamon/core/parallel_generic.h>
#include <boost/math/special_functions/legendre.hpp>
#include <math.h>
#include <pthread.h>
#include <stdexcept>
#include <time.h>
#include <tuple>
#include <vector>

// header for the deprecated code below
//#include <gsl/gsl_sf_legendre.h>
// this one is deprecated (GSL is REALLY slow)
// inline double spherical_harmonic(int l, int m, double theta, double phi){
//    int abs_m = (m>=0) ? m : -m;
//    double result = gsl_sf_legendre_sphPlm (l, abs_m, cos(theta));
//    if     (m > 0){
//        result *= sqrt2 * (((abs_m+1) % 2)-((abs_m) % 2)) * cos(abs_m*phi);
//    }
//    else{if(m < 0){
//        result *= sqrt2 * (((abs_m+1) % 2)-((abs_m) % 2)) * sin(abs_m*phi);
//    }}
//    return result;
//}
inline double spherical_harmonic(int l, int m, double theta, double phi) {
  int abs_m = (m >= 0) ? m : -m;
  double result = boost::math::legendre_p(l, abs_m, cos(theta)) *
                  sqrt_two_lplus1_div4pi[l] * sqrt_factorial[l - abs_m] *
                  one_div_sqrt_factorial[l + abs_m];
  if (m > 0) {
    result *= sqrt2 * (((abs_m + 1) % 2) - ((abs_m) % 2)) * cos(abs_m * phi);
  } else {
    if (m < 0) {
      result *= sqrt2 * (((abs_m + 1) % 2) - ((abs_m) % 2)) * sin(abs_m * phi);
    }
  }
  return result;
}

inline int binomial(int n, int k) {
  return (factorial[n]) / (factorial[n - k] * factorial[k]);
}

inline void cartesian_to_spherical(double x, double y, double z, double &r,
                                   double &theta, double &phi) {
  r = sqrt(x * x + y * y + z * z);
  theta = acos(z / r);
  phi = atan2(y, x);
  return;
}

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

std::ostream &operator<<(std::ostream &out, double val[3]) {
  out << "(" << val[0] << "|" << val[1] << "|" << val[2] << ")";
  return out;
}
std::ostream &operator<<(std::ostream &out, unsigned int val[3]) {
  out << "(" << val[0] << "|" << val[1] << "|" << val[2] << ")";
  return out;
}

inline double X_mu_nu(double A[3], double B[3], double P[3], unsigned int a[3],
                      unsigned int b[3], double d_mu_nu, double eta_mu_nu, RadInt RI,
                      const AngInt &AI) {
  double absP;
  double thetaP;
  double phiP;
  cartesian_to_spherical(P[0], P[1], P[2], absP, thetaP, phiP);

  RI.Init(eta_mu_nu, absP);

  ////LMAXP1 is defined in angular_integral.h
  // unsigned int lambda_max = a[0]+a[1]+a[2]+b[0]+b[1]+b[2];
  // double sph_harm[LMAXP1*LMAXP1] = {0};

  // double* sph_harm_it = sph_harm;
  // for (unsigned int lambda=0; lambda<=lambda_max; ++lambda){
  //    for (int gamma=-((int)lambda); gamma<=((int)(lambda)); ++gamma, ++sph_harm_it){

  //        *sph_harm_it = spherical_harmonic(lambda, gamma, thetaP, phiP);
  //    }
  //}

  double result = 0.0;
  for (unsigned int ax = 0; ax <= a[0]; ++ax /*, a_minus1=-a_minus1*/) {
    int a_binom[3];
    a_binom[0] = binomial(a[0], ax);
    for (unsigned int ay = 0; ay <= a[1]; ++ay /*, a_minus1=-a_minus1*/) {
      a_binom[1] = binomial(a[1], ay) * a_binom[0];
      for (unsigned int az = 0; az <= a[2]; ++az /*, a_minus1=-a_minus1*/) {
        a_binom[2] = binomial(a[2], az) * a_binom[1];

        int a_minus1 = ((a[0] + a[1] + a[2] - ax - ay - az) % 2 == 0) ? 1 : -1;
        unsigned int sum_a = ax + ay + az;
        double Apow = /*d_mu_nu **/ a_binom[2] * a_minus1 * power(A[0], a[0] - ax) *
                      power(A[1], a[1] - ay) * power(A[2], a[2] - az);

        for (unsigned int bx = 0; bx <= b[0]; ++bx /*, b_minus1=-b_minus1*/) {
          int b_binom[3];
          b_binom[0] = binomial(b[0], bx);
          for (unsigned int by = 0; by <= b[1]; ++by /*, b_minus1=-b_minus1*/) {
            b_binom[1] = binomial(b[1], by) * b_binom[0];
            for (unsigned int bz = 0; bz <= b[2]; ++bz /*, b_minus1=-b_minus1*/) {
              b_binom[2] = binomial(b[2], bz) * b_binom[1];

              int b_minus1 = ((b[0] + b[1] + b[2] - bx - by - bz) % 2 == 0) ? 1 : -1;
              unsigned int sum_b = bx + by + bz;
              double Bpow = b_binom[2] * b_minus1 * power(B[0], b[0] - bx) *
                            power(B[1], b[1] - by) * power(B[2], b[2] - bz);

              double Q = 0.0;
              // sph_harm_it = sph_harm;
              for (unsigned int lambda = 0; lambda <= sum_a + sum_b; ++lambda) {
                double Omega = 0.0;
                for (int gamma = -((int)lambda); gamma <= (int)(lambda);
                     ++gamma /*, ++sph_harm_it*/) {
                  double integral = AI.GetInt(lambda, gamma, ax + bx, ay + by, az + bz);
                  if (integral != 0.0) {
                    Omega += spherical_harmonic(lambda, gamma, thetaP, phiP) * integral;
                    // Omega += *sph_harm_it * integral;
                  }
                }
                if (Omega != 0.0) {
                  Q += RI.GetRadInt(sum_a + sum_b + 1, lambda) * Omega;
                  // Q += GetRadInt(eta_mu_nu, absP, sum_a+sum_b+1, lambda) * Omega;
                }
              }
              result += Apow * Bpow * Q;
            }
          }
        }
      }
    }
  }
  return result * d_mu_nu * 4.0 * Pi;
}

void *_potentialThreadOrbitals(void *data) {

  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
  struct timespec req = {0 /*req.tv_sec*/, 1L /*req.tv_nsec*/};
  GPSubData<double, double, double, int, double, double, double, double> *dat;
  dat = static_cast<
      GPSubData<double, double, double, int, double, double, double, double> *>(data);

  double *grdpnts = dat->GetData<0>();  // gridpoints
  double *prm_cent = dat->GetData<1>(); // centers of primitives
  int *prm_ang = dat->GetData<2>();     // angular exponents of primitives
  double *prm_exp = dat->GetData<3>();  // exponents of primitives
  double *prm_coef = dat->GetData<4>(); // contraction coefficients of primitives
  double *P_mu_nu = dat->GetData<5>();  // first order density matrix P_mu_nu
  double *S_mu_nu = dat->GetData<6>();  // S_mu_nu
  double *pot = dat->GetDataOutput();

  const int nr_prim = dat->GetNr<1>() / 3; // num of primitive gaussian basis functions
  const int nr_pot = dat->GetNrOutput();
  const int progress = 25;
  const bool progress_reports = dat->GetProgressReports();
  int *progress_bar = dat->GetProgressBar();
  pthread_mutex_t *mut = dat->GetMutex();

  RadInt RI = RadInt();
  const AngInt AI = AngInt();

  double *gp = grdpnts;
  double *potential = pot;
  for (int pp = 0; pp < nr_pot;) {
    for (int prog = 0; prog < progress && pp < nr_pot; ++prog, ++pp) {
      long double sum = 0.0;
      double r[3];
      r[0] = *gp;
      ++gp;
      r[1] = *gp;
      ++gp;
      r[2] = *gp;
      ++gp;
      // mu
      double *pcent_mu = prm_cent;
      int *pang_mu = prm_ang;
      double *pexp_mu = prm_exp;
      double *pcoef_mu = prm_coef;
      for (int mu = 0; mu < nr_prim; ++mu) {
        double A[3];
        A[0] = *pcent_mu - r[0];
        ++pcent_mu;
        A[1] = *pcent_mu - r[1];
        ++pcent_mu;
        A[2] = *pcent_mu - r[2];
        ++pcent_mu;
        unsigned int a[3];
        a[0] = (unsigned int)*pang_mu;
        ++pang_mu;
        a[1] = (unsigned int)*pang_mu;
        ++pang_mu;
        a[2] = (unsigned int)*pang_mu;
        ++pang_mu;
        double xi_mu = *pexp_mu;
        ++pexp_mu;
        double d_mu = *pcoef_mu;
        ++pcoef_mu;
        sum += X_mu_nu(A,
                       A,
                       A,
                       a,
                       a,
                       P_mu_nu[mu * nr_prim + mu] * d_mu * d_mu,
                       2.0 * xi_mu,
                       RI,
                       AI);
        // nu
        double *pcent_nu = pcent_mu;
        int *pang_nu = pang_mu;
        double *pexp_nu = pexp_mu;
        double *pcoef_nu = pcoef_mu;
        // both, P_mu_nu and X_mu_nu are symmetric
        for (int nu = mu + 1; nu < nr_prim; ++nu) {
          double B[3];
          B[0] = *pcent_nu - r[0];
          ++pcent_nu;
          B[1] = *pcent_nu - r[1];
          ++pcent_nu;
          B[2] = *pcent_nu - r[2];
          ++pcent_nu;
          unsigned int b[3];
          b[0] = (unsigned int)*pang_nu;
          ++pang_nu;
          b[1] = (unsigned int)*pang_nu;
          ++pang_nu;
          b[2] = (unsigned int)*pang_nu;
          ++pang_nu;
          double xi_nu = *pexp_nu;
          ++pexp_nu;
          double d_nu = *pcoef_nu;
          ++pcoef_nu;
          double eta_mu_nu = xi_mu + xi_nu;
          double P[3];
          P[0] = (A[0] - B[0]);
          P[1] = (A[1] - B[1]);
          P[2] = (A[2] - B[2]);
          double C_mu_nu =
              exp(-(xi_mu * xi_nu * ((P[0] * P[0]) + (P[1] * P[1]) + (P[2] * P[2]))) /
                  eta_mu_nu);
          P[0] = (xi_mu * A[0] + xi_nu * B[0]) / eta_mu_nu;
          P[1] = (xi_mu * A[1] + xi_nu * B[1]) / eta_mu_nu;
          P[2] = (xi_mu * A[2] + xi_nu * B[2]) / eta_mu_nu;
          double d_mu_nu = d_mu * d_nu * C_mu_nu * 2.0 * P_mu_nu[mu * nr_prim + nu];
          // WARNING: Increasing the cutoff for S_mu_nu or d_mu_nu might introduce
          // serious errors!!! WARNING: The current values screen away only such
          // integrals, that can be assumed to WARNING: definitely be insignificant. I
          // also tried 1.0e-10 and 1.0e-14, but these WARNING: caused serious
          // deviations from the real results.
          if (fabs(S_mu_nu[mu * nr_prim + nu]) > 1.0e-20 /*0.0*/ &&
              fabs(d_mu_nu) > 1.0e-30 /*0.0*/) {
            sum += X_mu_nu(A, B, P, a, b, d_mu_nu, eta_mu_nu, RI, AI);
          }
        }
      }
      *potential = sum;
      ++potential;
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

void electrostatic_potential_orbitals(
    bool progress_reports, int num_gridpoints, std::vector<double> prim_centers,
    std::vector<double> prim_exponents, std::vector<double> prim_coefficients,
    std::vector<int> prim_angular, std::vector<double> potential_grid,
    std::vector<double> P_matrix, std::vector<double> S_matrix,
    std::vector<double> *potential) {
  // initialize everything
  PG globals;
  init_parallel_generic(&progress_reports, &globals);

  const int split_factor = 3;

  // reserve data structures and fill them with input
  tuple_of_vectors<double, double, int, double, double, double, double> input;
  input = std::make_tuple(potential_grid,
                          prim_centers,
                          prim_angular,
                          prim_exponents,
                          prim_coefficients,
                          P_matrix,
                          S_matrix);

  // fill class that holds data for each thread
  GPData<double, double, double, int, double, double, double, double> *data;
  try {
    data = new GPData<double, double, double, int, double, double, double, double>(
        progress_reports,
        globals.nr_threads,
        input,
        potential,
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
  do_parallel_generic<double, double, double, int, double, double, double, double>(
      _potentialThreadOrbitals, &globals, progress_reports, num_gridpoints, data);
  // transfer output data
  data->TransferOutput();
  // clean up
  delete data;
  finalize_parallel_generic(progress_reports, &globals);
  pg_global = NULL;
}
