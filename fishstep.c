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

#include <pthread.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>

#if SDL
#include <SDL.h>
#endif

#if TIFF_ENABLE
#include <tiffio.h>
#endif

#include "fishstep.h"

// number of posix threads
#define NUM_THREADS 4

// draw every Nth frame
#define FRAMESKIP 2

// graphics scale up factor
#define SCALE 1

// number of cells per dimension
#define NUM_CELLS 512

// test problems, both 0: self-gravitating rotating disk
// SHOCK_TUBE 1: Sod shock tube along x, gas only, no gravity, writes shocktube.dat
// GAS_CLOUD 1: polytropic self-gravitating gas cylinder in hydrostatic
// equilibrium, writes radial profiles to gascloud.dat
#ifndef SHOCK_TUBE
#define SHOCK_TUBE 0
#endif
#ifndef GAS_CLOUD
#define GAS_CLOUD 0
#endif

#if SHOCK_TUBE
// high pressure gas in NUM_CELLS/4 <= x < 3*NUM_CELLS/4, low pressure outside,
// periodic box so there are two mirrored Riemann problems
// particles on a regular lattice, ST_PPC_*_X by ST_PPC_*_Y per cell,
// equal particle masses so the density ratio is the particle count ratio
#define ST_PPC_L_X 8
#define ST_PPC_L_Y 4
#define ST_PPC_R_X 2
#define ST_PPC_R_Y 2

// left state sound speed and pressure ratio P_R/P_L, Sod: rho_R/rho_L = 1/8
#define ST_CS_L 10.0
#define ST_P_RATIO 0.1

#define NUM_DM_PARTICLES 0
#define NUM_GAS_PARTICLES ((NUM_CELLS*NUM_CELLS/2)*(ST_PPC_L_X*ST_PPC_L_Y + ST_PPC_R_X*ST_PPC_R_Y))
#define GAS_MASS_FRACTION 1.0
#define GRAVITY 0
#ifndef T_END
#define T_END 5.0
#endif
#elif GAS_CLOUD
// gamma = 2 (n = 1) polytropic cylinder in hydrostatic equilibrium for
// lap(phi) = 4 pi G rho in 2D: P = K rho^2 gives 2K lap(rho) = -4 pi G rho,
// so rho(r) = rho_c J0(k r) with k^2 = 2 pi G/K and zero pressure at the
// radius R = j_01/k, the particles hold total mass 1 (a truncated isothermal
// cylinder is not an equilibrium in 2D, it is unbound below the critical
// line mass 2 cs^2/G)
#define GC_J01 2.404825557695773
#define GC_J1_AT_J01 0.5191474972894669
#define GC_R 48.0
#define GC_K_WAVE (GC_J01/GC_R)
#define GC_POLY_K (2.0*M_PI*G_EFF/(GC_K_WAVE*GC_K_WAVE))
#define GC_RHO_C (GC_K_WAVE/(2.0*M_PI*GC_R*GC_J1_AT_J01))
#define GC_OUT_DT 1.0

#define ISOTHERMAL 0
#define GAMMA 2.0

#define NUM_DM_PARTICLES 0
#define NUM_GAS_PARTICLES 1000000
#define GAS_MASS_FRACTION 1.0
#define GRAVITY 1
#ifndef T_END
#define T_END 10.0
#endif
#else
// number of particles per species
#define NUM_DM_PARTICLES (4*512000)
#define NUM_GAS_PARTICLES (1*512000)

// gas fraction of the total mass, total mass is normalized to 1
#define GAS_MASS_FRACTION 0.16

// self-gravity on, run until the window is closed
#define GRAVITY 1
#define T_END 0.0
#endif

#define NUM_PARTICLES (NUM_DM_PARTICLES + NUM_GAS_PARTICLES)

// display both species in separate color channels instead of total density
#define DISPLAY_SPECIES 1

// minimum and maximum geometries
#define XMIN 0
#define XMAX NUM_CELLS

// number of dimensions
#define NUM_DIMS 2

// particle mass per species
#define DM_MASS ((NUM_DM_PARTICLES) > 0 ? (1.0 - GAS_MASS_FRACTION)/(double)(NUM_DM_PARTICLES) : 0.0)
#define GAS_MASS ((NUM_GAS_PARTICLES) > 0 ? GAS_MASS_FRACTION/(double)(NUM_GAS_PARTICLES) : 0.0)

// size of the per-field, per-thread mesh subfield buffer
#define SUB_FIELD_SIZE (NUM_GRID_FIELDS*NUM_THREADS*NUM_CELLS*NUM_CELLS)

// gas equation of state
// ISOTHERMAL=1: P = CS_ISO^2 rho, internal energy is not evolved
// ISOTHERMAL=0: ideal gas P = (GAMMA - 1) rho u with adiabatic du/dt
#ifndef ISOTHERMAL
#define ISOTHERMAL 0
#endif
#ifndef GAMMA
#define GAMMA (5.0/3.0)
#endif
#ifndef CS_ISO
#define CS_ISO 5.0
#endif

// initial gas sound speed in cells per time unit, sets initial u
#define CS_INIT 5.0
#define U_INIT (CS_INIT*CS_INIT/(GAMMA*(GAMMA - 1.0)))

// internal energy floor
#define U_MIN 1.0e-6

// gas density below which a cell is treated as empty
#define RHO_GAS_FLOOR (1.0e-3*GAS_MASS)

// pressure forces act through cell faces, face f between cells i and j pushes
// them apart with F_f = lambda_f (P_i + P_j)/2. A sparse cell next to a dense
// one at a gas edge would be kicked far beyond the resolved acceleration, so
// lambda_f = min(1, rho_min/(RHO_ACC_EDGE*(rho_i + rho_j)/2)) limits the face
// and the force it no longer applies to the sparse cell stays with the dense
// one. Faces into empty cells have lambda_f = 0. Any lambda_f conserves
// momentum and energy exactly, the work of each face is taken from the
// internal energy of the two cells it pushes
#define RHO_ACC_EDGE 0.5

// von Neumann-Richtmyer artificial viscosity, linear and quadratic terms,
// applied in converging flow only
#define Q_LIN 0.5
#define Q_QUAD 2.0

// FLIP/PIC blending: each step gas particle velocities relax by this fraction
// towards the mass weighted mesh velocity, 0 is pure FLIP (particles keep
// their own velocities), 1 is pure PIC; damps sub-cell particle noise
#ifndef PIC_BLEND
#define PIC_BLEND 0.15
#endif

// return the kinetic energy removed by PIC blending to internal energy
#ifndef PIC_BLEND_HEAT
#define PIC_BLEND_HEAT 1
#endif

// maximum time step
#define DT 0.0095

// target CFL number for the adaptive time step, gas signal speed*dt/cell
#define CFL_TARGET 0.4

// effective gravitational constant of the mesh solver: the density source is
// scaled by 4 pi^2 and the unnormalized inverse FFT by 1/NUM_CELLS instead of
// 1/NUM_CELLS^2, so the potential solves lap(phi) = 4 pi G_EFF rho with unit cell size
#define G_EFF (M_PI*(double)(NUM_CELLS))

struct thread_data thread_data_array[NUM_THREADS];

// species table, dark matter first, gas after it
struct species_info species[NUM_SPECIES] = {
  [SPECIES_DM]  = {0, NUM_DM_PARTICLES, DM_MASS},
  [SPECIES_GAS] = {NUM_DM_PARTICLES, NUM_GAS_PARTICLES, GAS_MASS},
};

/* structure for RGB color */
struct color{
  double r;
  double g;
  double b;
};

void colormap(double val, struct color *col){
  int nn;

  double x_map[6]={0.0, 0.2, 0.45, 0.7, 0.85, 1.0};
  double r_map[6]={0.0, 0.0, 0.0, 1.0, 1.0, 0.65};
  double g_map[6]={0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
  double b_map[6]={0.65, 1.0, 1.0, 0.0, 0.0, 0.0};

  /* crop value to fit the map */
  if(val>0.99)
    val=0.99;

  val=1.0-val;
  
  /* linearly interpolate the value from the colormap */
  for(nn=0;nn<(6-1);nn++){
    if(val>x_map[nn]&&val<x_map[nn+1]){
      col->r=r_map[nn]+(val-x_map[nn])*(r_map[nn+1]-r_map[nn])/(x_map[nn+1]-x_map[nn]);
      col->g=g_map[nn]+(val-x_map[nn])*(g_map[nn+1]-g_map[nn])/(x_map[nn+1]-x_map[nn]);
      col->b=b_map[nn]+(val-x_map[nn])*(b_map[nn+1]-b_map[nn])/(x_map[nn+1]-x_map[nn]);
      break;
    }
  }
}

#if TIFF_ENABLE
void writeframe(char* path, SDL_Surface *screen){
  TIFF *file;
  Uint8 *p;
  int ii;

  int width=screen->w;
  int height=screen->h;

  file=TIFFOpen(path,"w");

  if(file){
    TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32) width);
    TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32) height);
    TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(file, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS);
    TIFFSetField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(file, TIFFTAG_EXTRASAMPLES, 0);
    TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField(file, TIFFTAG_IMAGEDESCRIPTION, "");

    p = (Uint8 *)screen->pixels;
    for (ii = height - 1; ii >= 0; ii--) {
      if (TIFFWriteScanline(file, p, ii, 0) < 0) {
	TIFFClose(file);
	printf("Error writing TIFF file.\n");
	exit(1);
      }
      p += 3 * width * sizeof(Uint8);
    }
    TIFFClose(file);
  }
}
#endif

// half kick of gas internal energy for the gas slice of a thread
static void kick_internal_energy(double *u, double *du, double dt, int thread_id){
#if !ISOTHERMAL
  int lo, hi, gg;

  thread_slice(0, NUM_GAS_PARTICLES, thread_id, NUM_THREADS, &lo, &hi);

  for(gg=lo; gg<hi; gg++) {
    u[gg] += 0.5*dt*du[gg];
    if(u[gg] < U_MIN) {
      u[gg] = U_MIN;
    }
  }
#endif
}

// particle kick-drift thread 
void *integration_kick_drift_thread(void *threadarg){
  // loop variables
  int ii;

  // thread slice boundaries
  int lo;
  int hi;

  // thread variables
  double *r, *v, *a;
  double *u, *du;
  double dt;

  // thread data pointer
  struct thread_data *my_data;

  // thread enumeration variable
  int thread_id;

  // set up pointers
  my_data=(struct thread_data *) threadarg;

  r=my_data->r;
  v=my_data->v;
  a=my_data->a;
  dt=my_data->dt;
  u=my_data->u;
  du=my_data->du;

  thread_id=my_data->thread_id;

  // compute thread slice
  lo=thread_id*NUM_PARTICLES/NUM_THREADS;
  hi=(thread_id+1)*NUM_PARTICLES/NUM_THREADS;

  for(ii=lo; ii<hi; ii++) {
    v[ii*NUM_DIMS + 0] += 0.5*dt*a[ii*NUM_DIMS + 0];
    v[ii*NUM_DIMS + 1] += 0.5*dt*a[ii*NUM_DIMS + 1];
    r[ii*NUM_DIMS + 0] += dt*v[ii*NUM_DIMS + 0];
    r[ii*NUM_DIMS + 1] += dt*v[ii*NUM_DIMS + 1];

    // enforce periodic boundaries
    r[ii*NUM_DIMS + 0] = wrap_periodic(r[ii*NUM_DIMS + 0], (double)(NUM_CELLS));
    r[ii*NUM_DIMS + 1] = wrap_periodic(r[ii*NUM_DIMS + 1], (double)(NUM_CELLS));
  }

  // half kick gas internal energy
  kick_internal_energy(u, du, dt, thread_id);

  pthread_exit(NULL);
}

// particle kick thread 
void *integration_kick_thread(void *threadarg){
  // loop variables
  int ii;

  // thread slice boundaries
  int lo;
  int hi;

  // thread variables
  double *v, *a;
  double *u, *du;
  double dt;

  // thread data pointer
  struct thread_data *my_data;

  // thread enumeration variable
  int thread_id;

  // set up pointers
  my_data=(struct thread_data *) threadarg;

  v=my_data->v;
  a=my_data->a;
  dt=my_data->dt;
  u=my_data->u;
  du=my_data->du;

  thread_id=my_data->thread_id;

  // compute thread slice
  lo=thread_id*NUM_PARTICLES/NUM_THREADS;
  hi=(thread_id+1)*NUM_PARTICLES/NUM_THREADS;

  for(ii=lo; ii<hi; ii++) {
    v[ii*NUM_DIMS + 0] += 0.5*dt*a[ii*NUM_DIMS + 0];
    v[ii*NUM_DIMS + 1] += 0.5*dt*a[ii*NUM_DIMS + 1];
  }

  // half kick gas internal energy
  kick_internal_energy(u, du, dt, thread_id);

  pthread_exit(NULL);
}

void *zero_rho_field_thread(void *threadarg){
  // loop variables
  int nn;

  // thread slice boundaries
  int lo;
  int hi;

  // thread variables
  double *sub_rho;
  
  // thread data pointer
  struct thread_data *my_data;

  // thread enumeration variable
  int thread_id;

  // set up pointers
  my_data=(struct thread_data *) threadarg;

  sub_rho=my_data->rho;

  thread_id=my_data->thread_id;

  // compute thread slice
  lo=thread_id*(SUB_FIELD_SIZE/NUM_THREADS);
  hi=(thread_id+1)*(SUB_FIELD_SIZE/NUM_THREADS);

  for(nn=lo; nn<hi; nn++) {
    sub_rho[nn] = 0.0;
  }

  pthread_exit(NULL);
}

void *copy_rho_to_complex_thread(void *threadarg){
  // loop variables
  int nn;

  // thread slice boundaries
  int lo;
  int hi;

  // thread variables
  fftw_complex *rho_complex;
  double *rho;
  
  // thread data pointer
  struct thread_data *my_data;

  // thread enumeration variable
  int thread_id;

  // set up pointers
  my_data=(struct thread_data *) threadarg;

  rho=my_data->rho;
  rho_complex=my_data->rho_complex;

  thread_id=my_data->thread_id;

  // compute thread slice
  lo=thread_id*(NUM_CELLS*NUM_CELLS)/NUM_THREADS;
  hi=(thread_id+1)*(NUM_CELLS*NUM_CELLS)/NUM_THREADS;

  for(nn=lo; nn<hi; nn++) {
    rho_complex[nn][0] = (double)(4.0*M_PI*M_PI)*rho[nn];
    rho_complex[nn][1] = 0.0;
  }

  pthread_exit(NULL);
}

void *compute_rho_field_thread(void *threadarg){
  // loop variables
  int nn;

  // thread slice boundaries
  int lo;
  int hi;

  // thread variables
  int ii,jj,kk;
  int s;
  double *r, *v, *u, *sub_rho, *rho;
  double *e_gas, *mx_gas, *my_gas;
  double dr[2];
  double w[4];
  double mass;
  int c[4];

  // thread data pointer
  struct thread_data *my_data;

  // thread enumeration variable
  int thread_id;

  // set up pointers
  my_data=(struct thread_data *) threadarg;

  r=my_data->r;
  v=my_data->v;
  u=my_data->u;
  sub_rho=my_data->rho;

  thread_id=my_data->thread_id;

  // gas energy and momentum subfields of this thread
  e_gas=&sub_rho[(FIELD_GAS_ENERGY*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];
  mx_gas=&sub_rho[(FIELD_GAS_MOM_X*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];
  my_gas=&sub_rho[(FIELD_GAS_MOM_Y*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];

  // deposit each species into its own per-thread subfield
  for(s=0; s<NUM_SPECIES; s++) {
    rho=&sub_rho[(s*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];
    mass=species[s].mass;

    // compute thread slice of this species
    thread_slice(species[s].offset, species[s].count, thread_id, NUM_THREADS, &lo, &hi);

    for(nn=lo; nn<hi; nn++) {
      ii = floor(r[nn*NUM_DIMS + 0]);
      jj = floor(r[nn*NUM_DIMS + 1]);

      // cloud in cell field interpolation
      dr[0] = r[nn*NUM_DIMS + 0] - (double)(ii);
      dr[1] = r[nn*NUM_DIMS + 1] - (double)(jj);

      c[0] = jj*NUM_CELLS + ii;
      c[1] = jj*NUM_CELLS + modulo(ii+1, NUM_CELLS);
      c[2] = modulo(jj+1, NUM_CELLS)*NUM_CELLS + ii;
      c[3] = modulo(jj+1, NUM_CELLS)*NUM_CELLS + modulo(ii+1, NUM_CELLS);

      w[0] = (1.0 - dr[0])*(1.0 - dr[1])*mass;
      w[1] = (dr[0])*(1.0 - dr[1])*mass;
      w[2] = (1.0 - dr[0])*(dr[1])*mass;
      w[3] = (dr[0])*(dr[1])*mass;

      for(kk=0; kk<4; kk++) {
	rho[c[kk]] += w[kk];
      }

      // gas also deposits internal energy and momentum
      if(s == SPECIES_GAS) {
	for(kk=0; kk<4; kk++) {
	  e_gas[c[kk]] += w[kk]*u[nn - species[SPECIES_GAS].offset];
	  mx_gas[c[kk]] += w[kk]*v[nn*NUM_DIMS + 0];
	  my_gas[c[kk]] += w[kk]*v[nn*NUM_DIMS + 1];
	}
      }
    }
  }

  pthread_exit(NULL);
}

void *acceleration_interpolation_thread(void *threadarg){
  // loop variables
  int nn;

  // thread slice boundaries
  int lo;
  int hi;

  // thread variables
  int ii,jj;
  double *r, *a, *phi;
  
  double dr[2];
  double g_x, g_y;

  // thread data pointer
  struct thread_data *my_data;

  // thread enumeration variable
  int thread_id;

  // set up pointers
  my_data=(struct thread_data *) threadarg;

  r=my_data->r;
  a=my_data->a;
  phi=my_data->phi;

  thread_id=my_data->thread_id;

  // compute thread slice
  lo=thread_id*NUM_PARTICLES/NUM_THREADS;
  hi=(thread_id+1)*NUM_PARTICLES/NUM_THREADS;

  for(nn=lo; nn<hi; nn++) {
    // calculate grid positions
    ii = floor(r[nn*NUM_DIMS + 0]);
    jj = floor(r[nn*NUM_DIMS + 1]);
    
    // calculate offset in grid
    dr[0] = r[nn*NUM_DIMS + 0] - (double)(ii);
    dr[1] = r[nn*NUM_DIMS + 1] - (double)(jj);

    // calculate force grid interpolation for 2nd order differences
    g_x = (phi[jj*NUM_CELLS + modulo(ii-1, NUM_CELLS)] - phi[jj*NUM_CELLS + modulo(ii+1, NUM_CELLS)])/2.0 * (1.0 - dr[0])*(1.0 - dr[1]);
    g_x += (phi[jj*NUM_CELLS + ii] - phi[jj*NUM_CELLS + modulo(ii+2, NUM_CELLS)])/2.0 * (dr[0])*(1.0-dr[1]);
    g_x += (phi[modulo(jj+1, NUM_CELLS)*NUM_CELLS + modulo(ii-1, NUM_CELLS)] - phi[modulo(jj+1, NUM_CELLS)*NUM_CELLS + modulo(ii+1, NUM_CELLS)])/2.0 * (1.0 - dr[0])*(dr[1]);
    g_x += (phi[modulo(jj+1, NUM_CELLS)*NUM_CELLS + ii] - phi[modulo(jj+1, NUM_CELLS)*NUM_CELLS + modulo(ii+2, NUM_CELLS)])/2.0 * (dr[0])*(dr[1]);

    g_y = (phi[modulo(jj-1, NUM_CELLS)*NUM_CELLS + ii] - phi[modulo(jj+1, NUM_CELLS)*NUM_CELLS + ii])/2.0 * (1.0 - dr[0])*(1.0 - dr[1]);
    g_y += (phi[modulo(jj-1, NUM_CELLS)*NUM_CELLS + modulo(ii+1, NUM_CELLS)] - phi[modulo(jj+1, NUM_CELLS)*NUM_CELLS + modulo(ii+1, NUM_CELLS)])/2.0 * (dr[0])*(1.0-dr[1]);
    g_y += (phi[jj*NUM_CELLS + ii] - phi[modulo(jj+2, NUM_CELLS)*NUM_CELLS + ii])/2.0 * (1.0 - dr[0])*(dr[1]);
    g_y += (phi[jj*NUM_CELLS + modulo(ii+1, NUM_CELLS)] - phi[modulo(jj+2, NUM_CELLS)*NUM_CELLS + modulo(ii+1, NUM_CELLS)])/2.0 * (dr[0])*(dr[1]);
    
    a[nn*NUM_DIMS + 0] = g_x;
    a[nn*NUM_DIMS + 1] = g_y;
  }

  pthread_exit(NULL);
}

// pressure force limiter of the face between cells with gas densities rho_i
// and rho_j, symmetric in its arguments so both cells see the same value
static inline double face_lambda(double rho_i, double rho_j){
  if(rho_i <= RHO_GAS_FLOOR || rho_j <= RHO_GAS_FLOOR) {
    return 0.0;
  }
  return fmin(1.0, fmin(rho_i, rho_j)/(RHO_ACC_EDGE*0.5*(rho_i + rho_j)));
}

// compute gas pressure including artificial viscosity, the face pressure
// acceleration and the matching specific internal energy rate on the mesh,
// returns maximum sound speed and velocity and the residual of the discrete
// work balance, which is zero up to roundoff
void compute_gas_pressure(const double *grid_fields, double *press, double *du_mesh, double *vel, double *acc, double *cs_max, double *v_max, double *p_res){
  int ii, jj, nn;
  int nxp, nxm, nyp, nym;
  double rho, cs2, divv, q, vabs;
  double l_xm, l_xp, l_ym, l_yp;

  const double *rho_gas = &grid_fields[FIELD_RHO_GAS*NUM_CELLS*NUM_CELLS];
  const double *e_gas = &grid_fields[FIELD_GAS_ENERGY*NUM_CELLS*NUM_CELLS];
  const double *mx_gas = &grid_fields[FIELD_GAS_MOM_X*NUM_CELLS*NUM_CELLS];
  const double *my_gas = &grid_fields[FIELD_GAS_MOM_Y*NUM_CELLS*NUM_CELLS];
  double *vx = &vel[0];
  double *vy = &vel[NUM_CELLS*NUM_CELLS];
  double *ax = &acc[0];
  double *ay = &acc[NUM_CELLS*NUM_CELLS];

  *cs_max = 0.0;
  *v_max = 0.0;

  // mass weighted mesh velocity and thermal pressure
  for(nn=0; nn<NUM_CELLS*NUM_CELLS; nn++) {
    rho = rho_gas[nn];
    if(rho > 0.0) {
      vx[nn] = mx_gas[nn]/rho;
      vy[nn] = my_gas[nn]/rho;
    }
    else {
      vx[nn] = 0.0;
      vy[nn] = 0.0;
    }

    if(rho > RHO_GAS_FLOOR) {
#if ISOTHERMAL
      press[nn] = CS_ISO*CS_ISO*rho;
      (void)e_gas;
#else
      press[nn] = (GAMMA - 1.0)*e_gas[nn];
#endif
    }
    else {
      press[nn] = 0.0;
    }
  }

  // velocity divergence and artificial viscosity
  for(jj=0; jj<NUM_CELLS; jj++) {
    for(ii=0; ii<NUM_CELLS; ii++) {
      nn = jj*NUM_CELLS + ii;
      rho = rho_gas[nn];

      if(rho <= RHO_GAS_FLOOR) {
	continue;
      }

      // central difference divergence, empty neighbors take the velocity of
      // this cell so the gas edge does not show up as spurious compression
      nxp = jj*NUM_CELLS + modulo(ii+1, NUM_CELLS);
      nxm = jj*NUM_CELLS + modulo(ii-1, NUM_CELLS);
      nyp = modulo(jj+1, NUM_CELLS)*NUM_CELLS + ii;
      nym = modulo(jj-1, NUM_CELLS)*NUM_CELLS + ii;

      divv = 0.5*(((rho_gas[nxp] > RHO_GAS_FLOOR) ? vx[nxp] : vx[nn]) -
		  ((rho_gas[nxm] > RHO_GAS_FLOOR) ? vx[nxm] : vx[nn])) +
	0.5*(((rho_gas[nyp] > RHO_GAS_FLOOR) ? vy[nyp] : vy[nn]) -
	     ((rho_gas[nym] > RHO_GAS_FLOOR) ? vy[nym] : vy[nn]));

#if ISOTHERMAL
      cs2 = CS_ISO*CS_ISO;
#else
      cs2 = GAMMA*press[nn]/rho;
#endif

      if(cs2 > *cs_max) {
	*cs_max = cs2;
      }
      vabs = vx[nn]*vx[nn] + vy[nn]*vy[nn];
      if(vabs > *v_max) {
	*v_max = vabs;
      }

      // viscous pressure in converging flow, cell size is 1
      q = 0.0;
      if(divv < 0.0) {
	q = rho*(Q_LIN*sqrt(cs2)*(-divv) + Q_QUAD*divv*divv);
      }

      press[nn] += q;
    }
  }

  // face pressure forces, acceleration f_c/rho_c and specific internal energy
  // rate at the mesh nodes. With lambda = 1 this is the central difference
  // -grad(P)/rho and -(P/rho) div(v). Each face exchanges equal and opposite
  // momentum between its cells, so sum_c rho_c a_c = 0, and the face work
  // F_f (V_j - V_i) gained as kinetic energy is removed from the cells in
  // proportion to their pressure, so sum_c (mom_c.a_c + rho_c du_c) = 0
  for(jj=0; jj<NUM_CELLS; jj++) {
    for(ii=0; ii<NUM_CELLS; ii++) {
      nn = jj*NUM_CELLS + ii;
      rho = rho_gas[nn];

      if(rho <= RHO_GAS_FLOOR) {
	ax[nn] = 0.0;
	ay[nn] = 0.0;
	du_mesh[nn] = 0.0;
	continue;
      }

      nxm = jj*NUM_CELLS + modulo(ii-1, NUM_CELLS);
      nxp = jj*NUM_CELLS + modulo(ii+1, NUM_CELLS);
      nym = modulo(jj-1, NUM_CELLS)*NUM_CELLS + ii;
      nyp = modulo(jj+1, NUM_CELLS)*NUM_CELLS + ii;

      l_xm = face_lambda(rho_gas[nxm], rho);
      l_xp = face_lambda(rho, rho_gas[nxp]);
      l_ym = face_lambda(rho_gas[nym], rho);
      l_yp = face_lambda(rho, rho_gas[nyp]);

      ax[nn] = 0.5*(l_xm*(press[nxm] + press[nn]) - l_xp*(press[nn] + press[nxp]))/rho;
      ay[nn] = 0.5*(l_ym*(press[nym] + press[nn]) - l_yp*(press[nn] + press[nyp]))/rho;

      du_mesh[nn] = -0.5*press[nn]*(l_xp*(vx[nxp] - vx[nn]) + l_xm*(vx[nn] - vx[nxm]) +
				    l_yp*(vy[nyp] - vy[nn]) + l_ym*(vy[nn] - vy[nym]))/rho;
    }
  }

  // residual of the discrete work balance: the pressure power on particles
  // sum_p m v.a = sum_c mom.a_c should cancel the heating sum_c rho du_c
  *p_res = 0.0;
  for(nn=0; nn<NUM_CELLS*NUM_CELLS; nn++) {
    *p_res += mx_gas[nn]*ax[nn] + my_gas[nn]*ay[nn] + rho_gas[nn]*du_mesh[nn];
  }

  *cs_max = sqrt(*cs_max);
  *v_max = sqrt(*v_max);
}

// add pressure acceleration to gas particles, compute du/dt and blend
// particle velocities towards the mesh velocity, all mesh quantities are
// interpolated with the CIC weights used for the deposit so that the
// pressure force and the blending conserve total momentum
void *hydro_interpolation_thread(void *threadarg){
  int nn, lo, hi, ii, jj, kk, gg;
  double *r, *v, *a, *u, *du;
  const double *ax, *ay, *vx, *vy, *du_mesh;
  double dr[2];
  double w[4];
  int c[4];
  double vpx, vpy, wx, wy, var, v2, h;
  double flip_loss = 0.0;

  struct thread_data *my_data;
  int thread_id;

  my_data=(struct thread_data *) threadarg;

  r=my_data->r;
  v=my_data->v;
  a=my_data->a;
  u=my_data->u;
  du=my_data->du;
  ax=&my_data->acc_mesh[0];
  ay=&my_data->acc_mesh[NUM_CELLS*NUM_CELLS];
  vx=&my_data->vel_mesh[0];
  vy=&my_data->vel_mesh[NUM_CELLS*NUM_CELLS];
  du_mesh=my_data->du_mesh;

  thread_id=my_data->thread_id;

  thread_slice(species[SPECIES_GAS].offset, species[SPECIES_GAS].count, thread_id, NUM_THREADS, &lo, &hi);

  for(nn=lo; nn<hi; nn++) {
    ii = floor(r[nn*NUM_DIMS + 0]);
    jj = floor(r[nn*NUM_DIMS + 1]);

    dr[0] = r[nn*NUM_DIMS + 0] - (double)(ii);
    dr[1] = r[nn*NUM_DIMS + 1] - (double)(jj);

    gg = nn - species[SPECIES_GAS].offset;

    // a = f/rho from the face pressure forces, divided by density on the mesh
    a[nn*NUM_DIMS + 0] += cic_interp(ax, ii, jj, dr, NUM_CELLS);
    a[nn*NUM_DIMS + 1] += cic_interp(ay, ii, jj, dr, NUM_CELLS);

    // du/dt from the face work, divided by density on the mesh
    du[gg] = cic_interp(du_mesh, ii, jj, dr, NUM_CELLS);

    if(PIC_BLEND > 0.0) {
      c[0] = jj*NUM_CELLS + ii;
      c[1] = jj*NUM_CELLS + modulo(ii+1, NUM_CELLS);
      c[2] = modulo(jj+1, NUM_CELLS)*NUM_CELLS + ii;
      c[3] = modulo(jj+1, NUM_CELLS)*NUM_CELLS + modulo(ii+1, NUM_CELLS);

      w[0] = (1.0 - dr[0])*(1.0 - dr[1]);
      w[1] = (dr[0])*(1.0 - dr[1]);
      w[2] = (1.0 - dr[0])*(dr[1]);
      w[3] = (dr[0])*(dr[1]);

      // mass weighted mesh velocity V at the particle and its CIC weighted
      // variance over the four cells
      vpx = 0.0;
      vpy = 0.0;
      for(kk=0; kk<4; kk++) {
	vpx += w[kk]*vx[c[kk]];
	vpy += w[kk]*vy[c[kk]];
      }
      var = 0.0;
      for(kk=0; kk<4; kk++) {
	var += w[kk]*((vx[c[kk]] - vpx)*(vx[c[kk]] - vpx) + (vy[c[kk]] - vpy)*(vy[c[kk]] - vpy));
      }

      // velocity relative to the mesh velocity, blended by a = PIC_BLEND
      wx = v[nn*NUM_DIMS + 0] - vpx;
      wy = v[nn*NUM_DIMS + 1] - vpy;

      v2 = v[nn*NUM_DIMS + 0]*v[nn*NUM_DIMS + 0] + v[nn*NUM_DIMS + 1]*v[nn*NUM_DIMS + 1];
      v[nn*NUM_DIMS + 0] -= PIC_BLEND*wx;
      v[nn*NUM_DIMS + 1] -= PIC_BLEND*wy;
      flip_loss += 0.5*(v2 - v[nn*NUM_DIMS + 0]*v[nn*NUM_DIMS + 0] - v[nn*NUM_DIMS + 1]*v[nn*NUM_DIMS + 1]);

      // the removed kinetic energy becomes heat. Per particle the blend removes
      // a(2-a)/2 |w|^2 + a V.w, summed over particles the second term is
      // a (sum_c rho_c |V_c|^2 - sum_p m |V_p|^2) because V_c = mom_c/rho_c,
      // which equals a sum_p m var_p. Heating each particle by
      // a(2-a)/2 |w|^2 + a var_p returns exactly the kinetic energy removed
      // with a local, Galilean invariant and non-negative heat
      if(PIC_BLEND_HEAT && !ISOTHERMAL) {
	h = 0.5*PIC_BLEND*(2.0 - PIC_BLEND)*(wx*wx + wy*wy) + PIC_BLEND*var;
	u[gg] += h;
	flip_loss -= h;
      }
    }
  }
  my_data->flip_loss = species[SPECIES_GAS].mass*flip_loss;

  pthread_exit(NULL);
}

// launch NUM_THREADS threads of fn on thread_data_array and join them
void run_threads(void *(*fn)(void *)){
  int nn, thread_rc;
  pthread_t threads[NUM_THREADS];
  pthread_attr_t attr;
  void *thread_status;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for(nn=0; nn<NUM_THREADS; nn++){
    thread_data_array[nn].thread_id=nn;
    thread_rc=pthread_create(&threads[nn], &attr, fn, (void *) &thread_data_array[nn]);
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
  }

  pthread_attr_destroy(&attr);
  for(nn=0; nn<NUM_THREADS; nn++){
    thread_rc=pthread_join(threads[nn], &thread_status);
    if(thread_rc){
      printf("ERROR: pthread_join() returned %d.\n", thread_rc);
      exit(-1);
    }
  }
}

#if GAS_CLOUD
// polytropic cylinder on a sunflower pattern: radii from the inverse of the
// enclosed mass fraction m(r) = r J1(k r)/(R J1(k R)) by bisection, golden
// angle azimuths, u = P/((gamma - 1) rho) = K rho for gamma = 2
void init_gas_cloud(double *r, double *v, double *a, double *u, double *du){
  int nn, it;
  double q, lo, hi, mid, th;
  double golden = M_PI*(3.0 - sqrt(5.0));

  for(nn=0; nn<NUM_GAS_PARTICLES; nn++) {
    q = ((double)(nn) + 0.5)/(double)(NUM_GAS_PARTICLES);
    lo = 0.0;
    hi = GC_R;
    for(it=0; it<60; it++) {
      mid = 0.5*(lo + hi);
      if(mid*j1(GC_K_WAVE*mid)/(GC_R*GC_J1_AT_J01) < q) {
	lo = mid;
      }
      else {
	hi = mid;
      }
    }
    th = golden*(double)(nn);

    r[nn*NUM_DIMS + 0] = 0.5*NUM_CELLS + mid*cos(th);
    r[nn*NUM_DIMS + 1] = 0.5*NUM_CELLS + mid*sin(th);
    v[nn*NUM_DIMS + 0] = 0.0;
    v[nn*NUM_DIMS + 1] = 0.0;
    a[nn*NUM_DIMS + 0] = 0.0;
    a[nn*NUM_DIMS + 1] = 0.0;
    u[nn] = GC_POLY_K*GC_RHO_C*j0(GC_K_WAVE*mid)/(GAMMA - 1.0);
    if(u[nn] < U_MIN) {
      u[nn] = U_MIN;
    }
    du[nn] = 0.0;
  }
}

static int compare_double(const void *x, const void *y){
  double dx = *(const double *)x, dy = *(const double *)y;
  return (dx > dy) - (dx < dy);
}

// append a radial density profile about the gas center of mass, normalized
// to the initial central density, and Lagrangian radii to fp
// total kinetic, gas thermal and gravitational energy of all particles into
// e[0], e[1], e[2], the potential energy W = 1/2 sum_c rho_c phi_c uses the
// mesh potential of the last solve
void compute_energy(const double *v, const double *u, const double *rho, const double *phi, double *e){
  int nn, s;
  double m;

  e[0] = 0.0;
  e[1] = 0.0;
  e[2] = 0.0;
  for(s=0; s<NUM_SPECIES; s++) {
    m = species[s].mass;
    for(nn=species[s].offset; nn<species[s].offset + species[s].count; nn++) {
      e[0] += 0.5*m*(v[nn*NUM_DIMS + 0]*v[nn*NUM_DIMS + 0] + v[nn*NUM_DIMS + 1]*v[nn*NUM_DIMS + 1]);
    }
  }
  for(nn=0; nn<NUM_GAS_PARTICLES; nn++) {
    e[1] += species[SPECIES_GAS].mass*u[nn];
  }
#if GRAVITY
  for(nn=0; nn<NUM_CELLS*NUM_CELLS; nn++) {
    e[2] += 0.5*rho[nn]*phi[nn];
  }
#else
  (void)rho;
  (void)phi;
#endif
}

// write the gas cloud radial profile block at time t, energy holds K, U, W,
// the cumulative FLIP kinetic energy loss not returned as heat and the
// cumulative residual of the discrete pressure work balance, or is NULL
// before the first mesh solve
void write_gas_cloud_profile(FILE *fp, const double *r, double t, const double *energy){
  int nn, bb;
  double cx = 0.0, cy = 0.0, dx, dy;
  double lagr_frac[5] = {0.1, 0.25, 0.5, 0.75, 0.9};
  double lagr[5];
  static double rad[NUM_GAS_PARTICLES];
  // radial bins of half a cell out to twice the truncation radius
  enum { NUM_BINS = (int)(4.0*GC_R) };
  double hist[NUM_BINS];

  for(nn=0; nn<NUM_GAS_PARTICLES; nn++) {
    cx += r[nn*NUM_DIMS + 0];
    cy += r[nn*NUM_DIMS + 1];
  }
  cx /= NUM_GAS_PARTICLES;
  cy /= NUM_GAS_PARTICLES;

  for(bb=0; bb<NUM_BINS; bb++) {
    hist[bb] = 0.0;
  }
  for(nn=0; nn<NUM_GAS_PARTICLES; nn++) {
    dx = r[nn*NUM_DIMS + 0] - cx;
    dy = r[nn*NUM_DIMS + 1] - cy;
    rad[nn] = sqrt(dx*dx + dy*dy);
    bb = (int)(2.0*rad[nn]);
    if(bb < NUM_BINS) {
      hist[bb] += GAS_MASS;
    }
  }

  qsort(rad, NUM_GAS_PARTICLES, sizeof(double), compare_double);
  for(bb=0; bb<5; bb++) {
    lagr[bb] = rad[(int)(lagr_frac[bb]*NUM_GAS_PARTICLES)];
  }

  fprintf(fp, "# t %.6f cm_x %.6f cm_y %.6f r10 %.4f r25 %.4f r50 %.4f r75 %.4f r90 %.4f R %.4f",
	  t, cx, cy, lagr[0], lagr[1], lagr[2], lagr[3], lagr[4], GC_R);
  if(energy) {
    fprintf(fp, " E_kin %.10e E_th %.10e E_pot %.10e flip_loss %.4e work_res %.4e",
	    energy[0], energy[1], energy[2], energy[3], energy[4]);
  }
  fprintf(fp, "\n");
  for(bb=0; bb<NUM_BINS; bb++) {
    // annulus area pi((r+0.5)^2 - r^2) with r = bb/2
    fprintf(fp, "%.4f %.8g\n", 0.5*bb + 0.25,
	    hist[bb]/(M_PI*0.25*(double)(2*bb + 1))/GC_RHO_C);
  }
  fprintf(fp, "\n\n");
  fflush(fp);

  printf("cloud t: %.3f center: (%.4f, %.4f) lagrangian radii r10 %.3f r25 %.3f r50 %.3f r75 %.3f r90 %.3f\n",
	 t, cx, cy, lagr[0], lagr[1], lagr[2], lagr[3], lagr[4]);
}
#endif

#if SHOCK_TUBE
// gas density, pressure and internal energy of the left (high pressure) state
#define ST_RHO_L (GAS_MASS*ST_PPC_L_X*ST_PPC_L_Y)
#define ST_RHO_R (GAS_MASS*ST_PPC_R_X*ST_PPC_R_Y)
#define ST_U_L (ST_CS_L*ST_CS_L/(GAMMA*(GAMMA - 1.0)))
#define ST_U_R (ST_U_L*ST_P_RATIO*ST_RHO_L/ST_RHO_R)

// Sod shock tube initial condition on a regular particle lattice
void init_shock_tube(double *r, double *v, double *a, double *u, double *du){
  int ii, jj, px, py, nx, ny;
  int nn = 0;
  int left;

  for(jj=0; jj<NUM_CELLS; jj++) {
    for(ii=0; ii<NUM_CELLS; ii++) {
      left = (ii >= NUM_CELLS/4 && ii < 3*NUM_CELLS/4);
      nx = left ? ST_PPC_L_X : ST_PPC_R_X;
      ny = left ? ST_PPC_L_Y : ST_PPC_R_Y;

      for(py=0; py<ny; py++) {
	for(px=0; px<nx; px++) {
	  r[nn*NUM_DIMS + 0] = ii + (px + 0.5)/nx;
	  r[nn*NUM_DIMS + 1] = jj + (py + 0.5)/ny;
	  v[nn*NUM_DIMS + 0] = 0.0;
	  v[nn*NUM_DIMS + 1] = 0.0;
	  a[nn*NUM_DIMS + 0] = 0.0;
	  a[nn*NUM_DIMS + 1] = 0.0;
	  u[nn] = left ? ST_U_L : ST_U_R;
	  du[nn] = 0.0;
	  nn++;
	}
      }
    }
  }

  if(nn != NUM_GAS_PARTICLES) {
    printf("Shock tube lattice placed %d particles, expected %d.\n", nn, NUM_GAS_PARTICLES);
    exit(1);
  }
}

// write y averaged gas profiles along x, normalized to the right state
// density and pressure, columns: x rho v_x P u
void write_shock_tube_profile(const char *path, const double *grid_fields, double t){
  int ii, jj, nn;
  double rho, mx, e;
  FILE *fp;

  const double *rho_gas = &grid_fields[FIELD_RHO_GAS*NUM_CELLS*NUM_CELLS];
  const double *e_gas = &grid_fields[FIELD_GAS_ENERGY*NUM_CELLS*NUM_CELLS];
  const double *mx_gas = &grid_fields[FIELD_GAS_MOM_X*NUM_CELLS*NUM_CELLS];
  double p_r = (GAMMA - 1.0)*ST_RHO_R*ST_U_R;

  fp = fopen(path, "w");
  if(!fp) {
    printf("Could not open %s for writing.\n", path);
    return;
  }

  // p_over_rho_R is the velocity unit squared of the normalized profile
  fprintf(fp, "# t %.10g gamma %.10g rho_L %.10g P_L %.10g x_L %d x_R %d p_over_rho_R %.10g\n",
	  t, GAMMA, ST_RHO_L/ST_RHO_R, (GAMMA - 1.0)*ST_RHO_L*ST_U_L/p_r, NUM_CELLS/4, 3*NUM_CELLS/4,
	  (GAMMA - 1.0)*ST_U_R);
  fprintf(fp, "# x rho v_x P u\n");

  for(ii=0; ii<NUM_CELLS; ii++) {
    rho = 0.0;
    mx = 0.0;
    e = 0.0;
    for(jj=0; jj<NUM_CELLS; jj++) {
      nn = jj*NUM_CELLS + ii;
      rho += rho_gas[nn];
      mx += mx_gas[nn];
      e += e_gas[nn];
    }
    // mesh nodes sit at integer coordinates
    fprintf(fp, "%d %.10g %.10g %.10g %.10g\n", ii,
	    rho/NUM_CELLS/ST_RHO_R, mx/rho, (GAMMA - 1.0)*e/NUM_CELLS/p_r, e/rho);
  }

  fclose(fp);
}
#endif

/* compute time difference in seconds and nanoseconds */
void timediff(struct timespec start, struct timespec end, struct timespec *out){
  /* compute time difference */
  if(end.tv_nsec<start.tv_nsec){
    out->tv_nsec=end.tv_nsec-start.tv_nsec+1000000000;
    out->tv_sec=end.tv_sec-start.tv_sec-1;
  }
  else{
    out->tv_nsec=end.tv_nsec-start.tv_nsec;
    out->tv_sec=end.tv_sec-start.tv_sec;
  }
}

int main(int argc, char *argv[])
{
  // loop variables
  int nn;
  int tt;
  double sim_time;
  int steps;
  int ii;
  int jj;
  
  // simulation variables
  double *r;
  double *v;
  double *a;

#if !SHOCK_TUBE && !GAS_CLOUD
  // initial condition variables
  double r2[NUM_DIMS];
  double dd;
  double vv;
#endif
  
  double *rho, *sub_rho, *rho_species;
  double *u, *du;
  double *press, *du_mesh, *gas_vel, *gas_acc;
  double cs_max, v_max, e_thermal;
  double p_res, e_flip_loss = 0.0, e_work_res = 0.0;
  double acc;
  int s;
  double *phi;
  double *gr_hat;

  double k_x, k_y, dd_k;
  
  fftw_complex *rho_hat, *green_rho_hat, *rho_complex, *phi_complex;
  fftw_plan rho_plan, phi_plan;
  
  // time step variable
  double dt=DT;

  struct color *col=0;

#if SDL
  // SDL variables
  SDL_Surface *screen;
  SDL_Event event;
  Uint8 *pixels;

  // SDL loop variables
  int ii2;
  int jj2;

  // gfx variables
  double d;
#if DISPLAY_SPECIES
  double d_dm, d_gas;
  double display_scale[NUM_SPECIES];
#endif
#endif

  // program state flag
  int done;

  // posix thread variables
  int thread_rc;
  pthread_t threads[NUM_THREADS];
  pthread_attr_t attr;
  void *thread_status;

  // time measurement variables
  struct timespec int_time;
  struct timespec time1, time2;

#if TIFF_ENABLE
  // tiff writer variables
  char filename[128];
  int tiff_frame;
#endif

  // init pseudorandom number generator
  srand(time(0));
  
#if SDL
  // open a SDL window
  SDL_Init(SDL_INIT_VIDEO);
  screen = SDL_SetVideoMode(SCALE*NUM_CELLS, SCALE*NUM_CELLS, 24, SDL_SWSURFACE);
  SDL_WM_SetCaption("Fishstep", "Fishstep");
#endif

  // allocate simulation vectors
  // displacement
  r=(double *)malloc(2*NUM_PARTICLES*sizeof(double));
  if(!r){
    printf("Out of memory: r not allocated.\n");
    exit(1);
  }
  // velocity
  v=(double *)malloc(2*NUM_PARTICLES*sizeof(double));
  if(!v){
    printf("Out of memory: v not allocated.\n");
    exit(1);
  }
  // acceleration
  a=(double *)malloc(2*NUM_PARTICLES*sizeof(double));
  if(!a){
    printf("Out of memory: a not allocated.\n");
    exit(1);
  }
  // density
  rho=(double *)malloc(NUM_CELLS*NUM_CELLS*sizeof(double));
  if(!rho){
    printf("Out of memory: rho not allocated.\n");
    exit(1);
  }
  sub_rho=(double *)malloc(SUB_FIELD_SIZE*sizeof(double));
  if(!sub_rho){
    printf("Out of memory: sub_rho not allocated.\n");
    exit(1);
  }
  // summed mesh fields, species densities first
  rho_species=(double *)malloc(NUM_GRID_FIELDS*NUM_CELLS*NUM_CELLS*sizeof(double));
  if(!rho_species){
    printf("Out of memory: rho_species not allocated.\n");
    exit(1);
  }
  // gas internal energy and its time derivative
  u=(double *)malloc(NUM_GAS_PARTICLES*sizeof(double));
  du=(double *)malloc(NUM_GAS_PARTICLES*sizeof(double));
  if(!u || !du){
    printf("Out of memory: u or du not allocated.\n");
    exit(1);
  }
  // gas pressure, du/dt and mesh velocity
  press=(double *)malloc(NUM_CELLS*NUM_CELLS*sizeof(double));
  du_mesh=(double *)malloc(NUM_CELLS*NUM_CELLS*sizeof(double));
  gas_vel=(double *)malloc(2*NUM_CELLS*NUM_CELLS*sizeof(double));
  gas_acc=(double *)malloc(2*NUM_CELLS*NUM_CELLS*sizeof(double));
  if(!press || !du_mesh || !gas_vel || !gas_acc){
    printf("Out of memory: gas pressure fields not allocated.\n");
    exit(1);
  }
  // potential
  phi=(double *)malloc(NUM_CELLS*NUM_CELLS*sizeof(double));
  if(!phi){
    printf("Out of memory: phi not allocated.\n");
    exit(1);
  }
  gr_hat=(double *)malloc(NUM_CELLS*NUM_CELLS*sizeof(double));
  if(!gr_hat){
    printf("Out of memory: gr_hat not allocated.\n");
    exit(1);
  }
  rho_complex = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NUM_CELLS*NUM_CELLS);
  if(!rho_complex){
    printf("Out of memory: rho_complex not allocated.\n");
    exit(1);
  }
  rho_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NUM_CELLS*NUM_CELLS);
  if(!rho_hat){
    printf("Out of memory: rho_hat not allocated.\n");
    exit(1);
  }
  phi_complex = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NUM_CELLS*NUM_CELLS);
  if(!phi_complex){
    printf("Out of memory: phi_complex not allocated.\n");
    exit(1);
  }
  green_rho_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NUM_CELLS*NUM_CELLS);
  if(!green_rho_hat){
    printf("Out of memory: green_rho_hat not allocated.\n");
    exit(1);
  }

  col=(struct color*)malloc(sizeof(struct color));
  if(!col){
    printf("Out of memory: col not allocated.\n");
    exit(1);
  }
  
  // generate initial condition
#if SHOCK_TUBE
  init_shock_tube(r, v, a, u, du);
#elif GAS_CLOUD
  init_gas_cloud(r, v, a, u, du);
#else
  // for now gas shares the collisionless rotating disk initial condition
  for(nn=0; nn<NUM_PARTICLES; nn++) {
    r[nn*NUM_DIMS + 0] = NUM_CELLS*(double)(rand())/RAND_MAX;
    r[nn*NUM_DIMS + 1] = NUM_CELLS*(double)(rand())/RAND_MAX;
    
    r2[0] = r[nn*NUM_DIMS + 0] - NUM_CELLS/2;
    r2[1] = r[nn*NUM_DIMS + 1] - NUM_CELLS/2;
    while(vector_norm(r2, 2) > 28.0){
      r[nn*NUM_DIMS + 0] = NUM_CELLS*(double)(rand())/RAND_MAX;
      r[nn*NUM_DIMS + 1] = NUM_CELLS*(double)(rand())/RAND_MAX;
      r2[0] = r[nn*NUM_DIMS + 0] - NUM_CELLS/2;
      r2[1] = r[nn*NUM_DIMS + 1] - NUM_CELLS/2;
    }

    dd = vector_norm(r2, 2);
    vv = sqrt(2.0*4.0*M_PI*M_PI*1.4825*dd);
    
    v[nn*NUM_DIMS + 0] = -r2[1]*vv/dd;
    v[nn*NUM_DIMS + 1] = r2[0]*vv/dd;

    //v[nn*NUM_DIMS + 0] = 0.0;
    //v[nn*NUM_DIMS + 1] = 0.0;

    a[nn*NUM_DIMS + 0] = 0.0;
    a[nn*NUM_DIMS + 1] = 0.0;
  }

  // gas internal energy
  for(nn=0; nn<NUM_GAS_PARTICLES; nn++) {
#if ISOTHERMAL
    u[nn] = CS_ISO*CS_ISO/(GAMMA - 1.0);
#else
    u[nn] = U_INIT;
#endif
    du[nn] = 0.0;
  }
#endif
  cs_max = 0.0;
  v_max = 0.0;
  
  // compute gravitational green's function in k-space for the 5 point
  // discrete laplacian, FFT bin ii holds wavenumber ii for ii <= N/2 and
  // ii - N above, the k = 0 mean density mode is dropped (periodic box)
  for(jj=0; jj<NUM_CELLS; jj++) {
    for(ii=0; ii<NUM_CELLS; ii++) {
      k_x = 2.0*M_PI*(double)(ii <= NUM_CELLS/2 ? ii : ii - NUM_CELLS)/(double)(NUM_CELLS);
      k_y = 2.0*M_PI*(double)(jj <= NUM_CELLS/2 ? jj : jj - NUM_CELLS)/(double)(NUM_CELLS);
      dd_k = sin(k_x/2.0)*sin(k_x/2.0) + sin(k_y/2.0)*sin(k_y/2.0);

      gr_hat[jj*NUM_CELLS + ii] = (dd_k > 0.0) ? -(1.0/4.0)/dd_k : 0.0;
    }
  }

#if FFTW3_THREADS  
  // initialize FFTW threads
  if(!fftw_init_threads()){
    printf("FFTW threads failed.\n");
    exit(1);
  }

  fftw_plan_with_nthreads(NUM_THREADS);
#endif
  
  // plan for fft2(rho)
  rho_plan = fftw_plan_dft_2d(NUM_CELLS, NUM_CELLS, rho_complex, rho_hat, FFTW_FORWARD, FFTW_ESTIMATE);
  if(!rho_plan){
    printf("FFTW rho_plan failed.\n");
    exit(1);
  }

  // plan for ifft2(gr_hat * rho_hat)
  phi_plan = fftw_plan_dft_2d(NUM_CELLS, NUM_CELLS, green_rho_hat, phi_complex, FFTW_BACKWARD, FFTW_ESTIMATE);
  if(!phi_plan){
    printf("FFTW phi_plan failed.\n");
    exit(1);
  }

  // start timer
  clock_gettime(CLOCK_MONOTONIC, &time1);

  // integrate time steps
  done=0;
  tt=0;
  sim_time=0.0;

#if GAS_CLOUD
  FILE *cloud_fp = fopen("gascloud.dat", "w");
  double next_out = GC_OUT_DT;
  double energy[5];
  if(!cloud_fp) {
    printf("Could not open gascloud.dat for writing.\n");
    exit(1);
  }
  printf("gas cloud: R %.2f cells, central sound speed %.4f cells/time, rho_c %.6e, crossing time R/cs_c %.4f\n",
	 GC_R, sqrt(GAMMA*GC_POLY_K*GC_RHO_C), GC_RHO_C, GC_R/sqrt(GAMMA*GC_POLY_K*GC_RHO_C));
  write_gas_cloud_profile(cloud_fp, r, 0.0, NULL);
#endif
  steps=0;
  
#if TIFF_ENABLE
  tiff_frame=0;
#endif

  while(!done){
#if SDL
    // check for SDL event
    if (SDL_PollEvent(&event)) {
      switch (event.type) {
        // close button clicked
      case SDL_QUIT:
	done=1;
	break;
      }
    }
#endif

    // main simulation loop
    // ********************

    // create kick-drift integration threads
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(nn=0;nn<NUM_THREADS;nn++){
      thread_data_array[nn].thread_id=nn;
      thread_data_array[nn].r=r;
      thread_data_array[nn].v=v;
      thread_data_array[nn].a=a;
      thread_data_array[nn].u=u;
      thread_data_array[nn].du=du;
      thread_data_array[nn].dt=dt;

      thread_rc=pthread_create(&threads[nn], &attr, integration_kick_drift_thread, (void *) &thread_data_array[nn]);
    
      if(thread_rc){
	printf("ERROR: pthread_create() returned %d.\n", thread_rc);
	exit(-1);
      }
    }
  
    /* join threads */
    pthread_attr_destroy(&attr);
    for(nn=0; nn < NUM_THREADS; nn++){
      thread_rc=pthread_join(threads[nn], &thread_status);
      
      if(thread_rc){
	printf("ERROR: pthread_join() returned %d.\n", thread_rc);
	exit(-1);
      }
    }

    // create zero thread rho field threads
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(nn=0;nn<NUM_THREADS;nn++){
      thread_data_array[nn].thread_id=nn;
      thread_data_array[nn].rho=sub_rho;

      thread_rc=pthread_create(&threads[nn], &attr, zero_rho_field_thread, (void *) &thread_data_array[nn]);
    
      if(thread_rc){
	printf("ERROR: pthread_create() returned %d.\n", thread_rc);
	exit(-1);
      }
    }
  
    /* join threads */
    pthread_attr_destroy(&attr);
    for(nn=0; nn < NUM_THREADS; nn++){
      thread_rc=pthread_join(threads[nn], &thread_status);
      
      if(thread_rc){
	printf("ERROR: pthread_join() returned %d.\n", thread_rc);
	exit(-1);
      }
    }

    // create rho field interpolation calculation threads
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(nn=0;nn<NUM_THREADS;nn++){
      thread_data_array[nn].thread_id=nn;
      thread_data_array[nn].r=r;
      thread_data_array[nn].v=v;
      thread_data_array[nn].u=u;
      thread_data_array[nn].rho=sub_rho;

      thread_rc=pthread_create(&threads[nn], &attr, compute_rho_field_thread, (void *) &thread_data_array[nn]);
    
      if(thread_rc){
	printf("ERROR: pthread_create() returned %d.\n", thread_rc);
	exit(-1);
      }
    }
  
    /* join threads */
    pthread_attr_destroy(&attr);
    for(nn=0;nn<NUM_THREADS;nn++){
      thread_rc=pthread_join(threads[nn], &thread_status);
      
      if(thread_rc){
	printf("ERROR: pthread_join() returned %d.\n", thread_rc);
	exit(-1);
      }
    }
    
    // sum thread subfields into mesh fields, species densities into total rho
    for(ii=0;ii<NUM_CELLS*NUM_CELLS;ii++){
      rho[ii] = 0.0;
      for(s=0;s<NUM_GRID_FIELDS;s++){
	acc = 0.0;
	for(nn=0;nn<NUM_THREADS;nn++){
	  acc += sub_rho[(s*NUM_THREADS + nn)*NUM_CELLS*NUM_CELLS + ii];
	}
	rho_species[s*NUM_CELLS*NUM_CELLS + ii] = acc;
	if(s < NUM_SPECIES){
	  rho[ii] += acc;
	}
      }
    }

#if GRAVITY
    // create rho to complex copy threads
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(nn=0;nn<NUM_THREADS;nn++){
      thread_data_array[nn].thread_id=nn;
      thread_data_array[nn].rho=rho;
      thread_data_array[nn].rho_complex=rho_complex;

      thread_rc=pthread_create(&threads[nn], &attr, copy_rho_to_complex_thread, (void *) &thread_data_array[nn]);
    
      if(thread_rc){
	printf("ERROR: pthread_create() returned %d.\n", thread_rc);
	exit(-1);
      }
    }
  
    /* join threads */
    pthread_attr_destroy(&attr);
    for(nn=0;nn<NUM_THREADS;nn++){
      thread_rc=pthread_join(threads[nn], &thread_status);
      
      if(thread_rc){
	printf("ERROR: pthread_join() returned %d.\n", thread_rc);
	exit(-1);
      }
    }
    
    // compute rho_hat
    fftw_execute(rho_plan);

    // multiply green's function and density in k-space
    for(jj=0; jj<NUM_CELLS; jj++) {
      for(ii=0; ii<NUM_CELLS; ii++) {
	green_rho_hat[jj*NUM_CELLS + ii][0] = gr_hat[jj*NUM_CELLS + ii] * rho_hat[jj*NUM_CELLS + ii][0];
	green_rho_hat[jj*NUM_CELLS + ii][1] = gr_hat[jj*NUM_CELLS + ii] * rho_hat[jj*NUM_CELLS + ii][1];
      }
    }
    
    // compute phi_complex
    fftw_execute(phi_plan);

    // compute real(phi_complex)
    for(ii=0; ii<NUM_CELLS*NUM_CELLS; ii++) {
      phi[ii] = 1.0/(NUM_CELLS) * phi_complex[ii][0];
    }

    // create acceleration interpolation threads
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(nn=0;nn<NUM_THREADS;nn++){
      thread_data_array[nn].thread_id=nn;
      thread_data_array[nn].r=r;
      thread_data_array[nn].a=a;
      thread_data_array[nn].phi=phi;
      thread_data_array[nn].dt=dt;

      thread_rc=pthread_create(&threads[nn], &attr, acceleration_interpolation_thread, (void *) &thread_data_array[nn]);
    
      if(thread_rc){
	printf("ERROR: pthread_create() returned %d.\n", thread_rc);
	exit(-1);
      }
    }
  
    /* join threads */
    pthread_attr_destroy(&attr);
    for(nn=0;nn<NUM_THREADS;nn++){
      thread_rc=pthread_join(threads[nn], &thread_status);
      
      if(thread_rc){
	printf("ERROR: pthread_join() returned %d.\n", thread_rc);
	exit(-1);
      }
    }
    
#else
    // no gravity, the hydro step adds pressure acceleration onto zero
    memset(a, 0, NUM_DIMS*NUM_PARTICLES*sizeof(double));
#endif

    // gas pressure on the mesh and pressure acceleration on gas particles
    compute_gas_pressure(rho_species, press, du_mesh, gas_vel, gas_acc, &cs_max, &v_max, &p_res);

    for(nn=0;nn<NUM_THREADS;nn++){
      thread_data_array[nn].r=r;
      thread_data_array[nn].v=v;
      thread_data_array[nn].a=a;
      thread_data_array[nn].u=u;
      thread_data_array[nn].du=du;
      thread_data_array[nn].grid_fields=rho_species;
      thread_data_array[nn].acc_mesh=gas_acc;
      thread_data_array[nn].vel_mesh=gas_vel;
      thread_data_array[nn].du_mesh=du_mesh;
    }
    run_threads(hydro_interpolation_thread);
    for(nn=0;nn<NUM_THREADS;nn++){
      e_flip_loss += thread_data_array[nn].flip_loss;
    }
    e_work_res += p_res*dt;

    // create kick integration threads
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(nn=0;nn<NUM_THREADS;nn++){
      thread_data_array[nn].thread_id=nn;
      thread_data_array[nn].r=r;
      thread_data_array[nn].v=v;
      thread_data_array[nn].a=a;
      thread_data_array[nn].u=u;
      thread_data_array[nn].du=du;
      thread_data_array[nn].dt=dt;

      thread_rc=pthread_create(&threads[nn], &attr, integration_kick_thread, (void *) &thread_data_array[nn]);
    
      if(thread_rc){
	printf("ERROR: pthread_create() returned %d.\n", thread_rc);
	exit(-1);
      }
    }
  
    /* join threads */
    pthread_attr_destroy(&attr);
    for(nn=0;nn<NUM_THREADS;nn++){
      thread_rc=pthread_join(threads[nn], &thread_status);
      
      if(thread_rc){
	printf("ERROR: pthread_join() returned %d.\n", thread_rc);
	exit(-1);
      }
    }

    // end main simulation loop
    // ************************
    
    // draw state
    if(steps == FRAMESKIP-1) {
    
      // timer stop
      clock_gettime(CLOCK_MONOTONIC, &time2);

#if SDL
      pixels=(Uint8 *)screen->pixels;
      SDL_LockSurface(screen);

#if DISPLAY_SPECIES
      // normalize each species by its mass fraction so a well mixed
      // species has the same brightness as the total density display
      for(s=0;s<NUM_SPECIES;s++){
	acc = species[s].mass*(double)(species[s].count);
	display_scale[s] = (acc > 0.0) ? NUM_PARTICLES*4.0/acc : 0.0;
      }
#endif
    
      for(jj=0; jj<NUM_CELLS; jj++){
	for(ii=0; ii<NUM_CELLS; ii++){
	  // get field value
	  d=NUM_PARTICLES*4.0*rho[jj*NUM_CELLS+ii];
#if DISPLAY_SPECIES
	  d_dm=display_scale[SPECIES_DM]*rho_species[SPECIES_DM*NUM_CELLS*NUM_CELLS + jj*NUM_CELLS+ii];
	  d_gas=display_scale[SPECIES_GAS]*rho_species[SPECIES_GAS*NUM_CELLS*NUM_CELLS + jj*NUM_CELLS+ii];
	  d_dm=fmin(fmax(d_dm, 0.0), 255.0);
	  d_gas=fmin(fmax(d_gas, 0.0), 255.0);
#endif
	  //d=-0.1*phi[jj*NUM_CELLS+ii];

	  // clip value
	  if(d>255.0) {
	    d=255.0;
	  }
	  else if(d<0.0) {
	    d=0.0;
	  }

	  //d=d/256.0;
	  
	  // update framebuffer
	  for(jj2=0; jj2<SCALE; jj2++){
	    for(ii2=0; ii2<SCALE; ii2++){
	      // display heatmap style field
	      //colormap(d, col);
	      //pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+0]=(Uint8)(col->r*255.0);
	      //pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+1]=(Uint8)(col->g*255.0);
	      //pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+2]=(Uint8)(col->b*255.0);

#if DISPLAY_SPECIES
	      // display gas and collisionless densities in separate channels
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+0]=(Uint8)(d_gas);
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+1]=(Uint8)(0.5*(d_gas+d_dm));
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+2]=(Uint8)(d_dm);
#else
	      // display density of particles
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+0]=(Uint8)(d);
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+1]=(Uint8)(d);
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+2]=(Uint8)(d);
#endif
	    }
	  }
	}
      }

      SDL_UnlockSurface(screen);
      SDL_Flip(screen);
#endif
    
      timediff(time1, time2, &int_time);

      // gas thermal energy
      e_thermal = 0.0;
      for(nn=0; nn<NUM_GAS_PARTICLES; nn++) {
	e_thermal += GAS_MASS*u[nn];
      }

      // print state, CFL is the gas signal speed times dt over cell size
      printf("tt: %d t: %.4f dt: %.3e steps/s: %d E_th: %.6e cs_max: %.3f v_max: %.3f CFL: %.3f\n", tt,
	     sim_time, dt,
	     (int)(steps/((double)(int_time.tv_sec)+(double)(int_time.tv_nsec)*1.0E-9)),
	     e_thermal, cs_max, v_max, (cs_max + v_max)*dt);

#if TIFF_ENABLE
      sprintf(filename, "/home/janne808/testrun/%08d.tif", tiff_frame++);
      writeframe(filename, screen);
#endif

      // timer start
      clock_gettime(CLOCK_MONOTONIC, &time1);

      // reset steps variable
      steps = -1;
    }
    // adaptive time step for the next step from the gas CFL condition,
    // changed only between full steps so both half kicks use the same dt
    sim_time += dt;
    dt = DT;
    if((cs_max + v_max)*DT > CFL_TARGET) {
      dt = CFL_TARGET/(cs_max + v_max);
    }

#if GAS_CLOUD
    // profile output, the step is shortened to land on output times
    if(sim_time >= next_out - 1.0e-12) {
      compute_energy(v, u, rho, phi, energy);
      energy[3] = e_flip_loss;
      energy[4] = e_work_res;
      printf("energy t: %.4f K: %.8e U: %.8e W: %.8e E: %.8e flip loss: %.3e work residual: %.3e\n",
	     sim_time, energy[0], energy[1], energy[2], energy[0] + energy[1] + energy[2], energy[3], energy[4]);
      write_gas_cloud_profile(cloud_fp, r, sim_time, energy);
      next_out += GC_OUT_DT;
    }
    if(sim_time + dt > next_out) {
      dt = next_out - sim_time;
    }
#endif

    // stop at T_END, last step is shortened to land on it
    if(T_END > 0.0) {
      if(sim_time >= T_END - 1.0e-12) {
#if SHOCK_TUBE
	write_shock_tube_profile("shocktube.dat", rho_species, sim_time);
#endif
	printf("reached T_END at tt: %d t: %.6f\n", tt, sim_time);
	done = 1;
      }
      else if(sim_time + dt > T_END) {
	dt = T_END - sim_time;
      }
    }

    // update time step variable
    tt++;
    steps++;
  }

#if GAS_CLOUD
  fclose(cloud_fp);
#endif

  // clean up FFTW
  fftw_destroy_plan(rho_plan);
#if FFTW3_THREADS  
  fftw_cleanup_threads();
#endif
  
#if SDL
  // clean up SDL
  SDL_Quit();
#endif

  return 0;
}

