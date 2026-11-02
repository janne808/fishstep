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
  SPECIES_DUST = 2, // pressureless dust superparticles coupled to the gas by drag
  NUM_SPECIES = 3
};

// mesh fields deposited from particles with cloud in cell weights
// the first NUM_SPECIES fields are the species mass densities
enum grid_field_id{
  FIELD_RHO_DM = SPECIES_DM,
  FIELD_RHO_GAS = SPECIES_GAS,
  FIELD_RHO_DUST = SPECIES_DUST,
  FIELD_GAS_ENERGY = NUM_SPECIES,  // gas internal energy density, sum of m*u
  FIELD_GAS_MOM_X,                 // gas momentum density, sum of m*v_x
  FIELD_GAS_MOM_Y,                 // gas momentum density, sum of m*v_y
  FIELD_GAS_AFF_X,                 // APIC affine momentum, sum of m*(B grad(w))_x
  FIELD_GAS_AFF_Y,                 // APIC affine momentum, sum of m*(B grad(w))_y
  FIELD_DUST_MOM_X,                // dust momentum density, sum of m*v_x
  FIELD_DUST_MOM_Y,                // dust momentum density, sum of m*v_y
  FIELD_DUST_AFF_X,                // dust APIC affine momentum, sum of m*(B grad(w))_x
  FIELD_DUST_AFF_Y,                // dust APIC affine momentum, sum of m*(B grad(w))_y
  NUM_GRID_FIELDS
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

  // gas specific internal energy and its time derivative, indexed by
  // gas particle number (particle index minus gas species offset)
  double *u;
  double *du;

  // APIC affine velocity matrix of each gas particle, B[4*g + 2*a + b] holds
  // B_ab = sum_c w_c V_c,a (x_c - x_p)_b, the particle velocity field around
  // x_p is v_p + B D^-1 (x - x_p) with D = sum_c w_c (x_c - x_p)(x_c - x_p)^T
  double *aff;
  // the same for dust particles, indexed by dust particle number
  double *daff;

  // nonzero for gas particles accreted by a sink (the star of the disk
  // problem): they move with it and deposit into the collisionless density
  unsigned char *sink;
  // the same for dust particles, indexed by dust particle number
  unsigned char *dsink;

  // dust drag on the mesh (see compute_dust_drag) and the per-thread
  // subfields that collect the kinetic energy change of the drag kick
  double *drag_mesh;
  double *sub_heat;

  // summed mesh fields, gas pressure acceleration, mass weighted gas
  // velocity including the APIC affine momentum (first two of the four
  // velocity components) and specific internal energy rate on the mesh
  double *grid_fields;
  double *acc_mesh;
  double *vel_mesh;
  double *du_mesh;

  // kinetic energy removed by FLIP blending minus the heat returned, per thread
  double flip_loss;

  // kinetic energy change of the drag kick and the heat returned, per thread
  double drag_ke;
  double drag_heat;

  double dt;
};

// vector norm
static inline double vector_norm(double *r, int dim) {
  int ii;
  double s=0;
  
  for(ii=0; ii<dim; ii++) {
    s += r[ii]*r[ii];
  }

  return sqrt(s);
}

// custom modulo function to patch C style % operator
static inline int modulo(int x, int mod){
  if(x >= 0 && x < mod){
    return x;
  }
  else{
    return ((x % mod) + mod) % mod;
  }
}

// wrap a coordinate into the periodic interval [0, L)
static inline double wrap_periodic(double x, double L){
  x = fmod(x, L);
  if(x < 0.0) {
    x += L;
  }
  // a tiny negative x rounds to exactly L after the shift
  if(x >= L) {
    x = 0.0;
  }
  return x;
}

// cloud in cell stencil of a particle at offset dr inside cell (ii, jj): the
// four node indices, weights w, weight gradients gw and node offsets
// d = x_c - x_p, nodes ordered (ii, jj), (ii+1, jj), (ii, jj+1), (ii+1, jj+1)
static inline void cic_stencil(int ii, int jj, const double *dr, int n,
			       int *c, double *w, double gw[4][2], double d[4][2]){
  int ip = modulo(ii+1, n);
  int jp = modulo(jj+1, n);

  c[0] = jj*n + ii;
  c[1] = jj*n + ip;
  c[2] = jp*n + ii;
  c[3] = jp*n + ip;

  w[0] = (1.0 - dr[0])*(1.0 - dr[1]);
  w[1] = (dr[0])*(1.0 - dr[1]);
  w[2] = (1.0 - dr[0])*(dr[1]);
  w[3] = (dr[0])*(dr[1]);

  gw[0][0] = -(1.0 - dr[1]); gw[0][1] = -(1.0 - dr[0]);
  gw[1][0] = (1.0 - dr[1]);  gw[1][1] = -dr[0];
  gw[2][0] = -dr[1];         gw[2][1] = (1.0 - dr[0]);
  gw[3][0] = dr[1];          gw[3][1] = dr[0];

  d[0][0] = -dr[0];       d[0][1] = -dr[1];
  d[1][0] = 1.0 - dr[0];  d[1][1] = -dr[1];
  d[2][0] = -dr[0];       d[2][1] = 1.0 - dr[1];
  d[3][0] = 1.0 - dr[0];  d[3][1] = 1.0 - dr[1];
}

// cloud in cell interpolation of a periodic field f to offset dr inside cell (ii, jj)
static inline double cic_interp(const double *f, int ii, int jj, const double *dr, int n){
  int ip = modulo(ii+1, n);
  int jp = modulo(jj+1, n);

  return f[jj*n + ii]*(1.0 - dr[0])*(1.0 - dr[1]) +
    f[jj*n + ip]*(dr[0])*(1.0 - dr[1]) +
    f[jp*n + ii]*(1.0 - dr[0])*(dr[1]) +
    f[jp*n + ip]*(dr[0])*(dr[1]);
}

#endif
