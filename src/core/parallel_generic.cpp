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
#include <FireDeamon/core/parallel_generic.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>

// PG stands for ParallelGlobals
PG::PG() {
  this->threads = NULL;
  this->nr_threads = 0;
  this->progress_bar = 0;
}

PG::~PG() {
  if (this->threads != NULL) {
    free(this->threads);
  }
}

PG *pg_global = NULL;

void signal_callback_handler(int signum) {
  fprintf(stderr, "\nCaught signal %d\n", signum);
  if (pg_global != NULL) {
    if (pg_global->threads != NULL) {
      for (int thread_count = 0; thread_count < pg_global->nr_threads; ++thread_count) {
        int rc = pthread_cancel(*(pg_global->threads + thread_count));
        if (rc) {
          fprintf(
              stderr,
              "Cancellation of thread number %d out of %d could not be initialized.\n",
              thread_count + 1,
              pg_global->nr_threads);
        }
      }
    }
  }
  exit(signum);
}

void init_parallel_generic(bool *progress_reports, PG *globals) {
  int rc;
  int num_threads = 1;
  {
    char const *tmp = getenv("OMP_NUM_THREADS");
    if (tmp != NULL) {
      num_threads = atoi(tmp);
    }
  }
  globals->nr_threads = num_threads;
  globals->threads = new pthread_t[num_threads];
  globals->progress_bar = 0;
  if (*progress_reports) {
    if (globals->nr_threads > 1) {
      printf("Starting multi threaded computation with %d threads.\n",
             globals->nr_threads);
    } else {
      printf("Starting single threaded computation.\n");
    }
    fflush(stdout);
    rc = pthread_mutex_init(&(globals->mutex), NULL);
    if (rc) {
      fprintf(stderr, "Failed to initialize mutex, will disable progress reports.\n");
      *progress_reports = false;
    }
  }
}

void finalize_parallel_generic(bool progress_reports, PG *globals) {
  if (progress_reports) {
    int rc = pthread_mutex_trylock(&(globals->mutex));
    if (rc) {
      fprintf(stderr, "The mutex could not be acquired, terminating abnormally.\n");
    } else {
      pthread_mutex_unlock(&(globals->mutex));
      pthread_mutex_destroy(&(globals->mutex));
    }
  }
}
