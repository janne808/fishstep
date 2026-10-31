/*
 *  (C) 2020 Janne Heikkarainen <janne808@radiofreerobotron.net>
 *
 *  All rights reserved.
 *
 *  This file is part of Fishstep N-body Simulator.
 *
 *  Fishstep N-body Simulator is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Fishstep N-body Simulator is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Fishstep N-body Simulator.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __FISHSTEPH__
#define __FISHSTEPH__

// particle species
// particles of each species are stored contiguously in the shared r, v, a
// arrays: species s occupies indices [offset, offset + count)
enum species_id{
  SPECIES_DM = 0,   // collisionless (dark matter / stars)
  SPECIES_GAS = 1,  // gas, to be coupled to PIC hydrodynamics
  NUM_SPECIES = 2
};

struct species_info{
  int offset;
  int count;
  double mass;
};

// split a contiguous particle range [offset, offset + count) into
// num_threads slices and return the slice of thread_id in [*lo, *hi)
static inline void thread_slice(int offset, int count, int thread_id, int num_threads, int *lo, int *hi){
  *lo = offset + (int)((long)thread_id*count/num_threads);
  *hi = offset + (int)((long)(thread_id + 1)*count/num_threads);
}

struct thread_data{
  int thread_id;

  double *r;
  double *v;
  double *a;

  double *phi;
  double *rho;
  fftw_complex *rho_complex;
  
  double dt;
};

// vector norm
inline double vector_norm(double *r, int dim) {
  int ii;
  double s=0;
  
  for(ii=0; ii<dim; ii++) {
    s += r[ii]*r[ii];
  }

  return sqrt(s);
}

// custom modulo function to patch C style % operator
inline int modulo(int x, int mod){
  if(x >= 0 && x < mod){
    return x;
  }
  else{
    return ((x % mod) + mod) % mod;
  }
}

#endif
