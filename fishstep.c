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

// test problem, 0: self-gravitating rotating disk
// 1: Sod shock tube along x, gas only, no gravity, writes shocktube.dat at T_END
#ifndef SHOCK_TUBE
#define SHOCK_TUBE 0
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
#define T_END 5.0
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
#define ISOTHERMAL 0
#define GAMMA (5.0/3.0)
#define CS_ISO 5.0

// initial gas sound speed in cells per time unit, sets initial u
#define CS_INIT 5.0
#define U_INIT (CS_INIT*CS_INIT/(GAMMA*(GAMMA - 1.0)))

// internal energy floor
#define U_MIN 1.0e-6

// gas density below which a cell is treated as empty
#define RHO_GAS_FLOOR (1.0e-3*GAS_MASS)

// von Neumann-Richtmyer artificial viscosity, linear and quadratic terms,
// applied in converging flow only
#define Q_LIN 0.5
#define Q_QUAD 2.0

// maximum time step
#define DT 0.0095

// target CFL number for the adaptive time step, gas signal speed*dt/cell
#define CFL_TARGET 0.4

#define G 4.0*M_PI*M_PI;

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

// compute gas pressure including artificial viscosity and specific internal
// energy rate -(P/rho) div(v) on the mesh, returns maximum sound speed and velocity
void compute_gas_pressure(const double *grid_fields, double *press, double *du_mesh, double *vel, double *cs_max, double *v_max){
  int ii, jj, nn;
  int nxp, nxm, nyp, nym;
  double rho, cs2, divv, q, vabs;

  const double *rho_gas = &grid_fields[FIELD_RHO_GAS*NUM_CELLS*NUM_CELLS];
  const double *e_gas = &grid_fields[FIELD_GAS_ENERGY*NUM_CELLS*NUM_CELLS];
  const double *mx_gas = &grid_fields[FIELD_GAS_MOM_X*NUM_CELLS*NUM_CELLS];
  const double *my_gas = &grid_fields[FIELD_GAS_MOM_Y*NUM_CELLS*NUM_CELLS];
  double *vx = &vel[0];
  double *vy = &vel[NUM_CELLS*NUM_CELLS];

  *cs_max = 0.0;
  *v_max = 0.0;

  // mesh velocity and thermal pressure
  for(nn=0; nn<NUM_CELLS*NUM_CELLS; nn++) {
    rho = rho_gas[nn];
    if(rho > RHO_GAS_FLOOR) {
      vx[nn] = mx_gas[nn]/rho;
      vy[nn] = my_gas[nn]/rho;
#if ISOTHERMAL
      press[nn] = CS_ISO*CS_ISO*rho;
#else
      press[nn] = (GAMMA - 1.0)*e_gas[nn];
#endif
    }
    else {
      vx[nn] = 0.0;
      vy[nn] = 0.0;
      press[nn] = 0.0;
    }
  }

  // velocity divergence, artificial viscosity and p div v
  for(jj=0; jj<NUM_CELLS; jj++) {
    for(ii=0; ii<NUM_CELLS; ii++) {
      nn = jj*NUM_CELLS + ii;
      rho = rho_gas[nn];

      if(rho <= RHO_GAS_FLOOR) {
	du_mesh[nn] = 0.0;
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

      // specific internal energy rate -(P/rho) div v
      du_mesh[nn] = -press[nn]*divv/rho;
    }
  }

  *cs_max = sqrt(*cs_max);
  *v_max = sqrt(*v_max);
}

// add pressure acceleration to gas particles and compute du/dt
void *hydro_interpolation_thread(void *threadarg){
  int nn, lo, hi, ii, jj;
  double *r, *a, *du;
  const double *rho_gas, *press, *du_mesh;
  double dr[2];
  double gx, gy, rho_p;

  struct thread_data *my_data;
  int thread_id;

  my_data=(struct thread_data *) threadarg;

  r=my_data->r;
  a=my_data->a;
  du=my_data->du;
  rho_gas=&my_data->grid_fields[FIELD_RHO_GAS*NUM_CELLS*NUM_CELLS];
  press=my_data->press;
  du_mesh=my_data->du_mesh;

  thread_id=my_data->thread_id;

  thread_slice(species[SPECIES_GAS].offset, species[SPECIES_GAS].count, thread_id, NUM_THREADS, &lo, &hi);

  for(nn=lo; nn<hi; nn++) {
    ii = floor(r[nn*NUM_DIMS + 0]);
    jj = floor(r[nn*NUM_DIMS + 1]);

    dr[0] = r[nn*NUM_DIMS + 0] - (double)(ii);
    dr[1] = r[nn*NUM_DIMS + 1] - (double)(jj);

    // gas density at particle, nonzero since the particle deposited itself
    rho_p = cic_interp(rho_gas, ii, jj, dr, NUM_CELLS);

    if(rho_p > 0.0) {
      // a = -grad(P)/rho
      cic_neg_gradient(press, ii, jj, dr, NUM_CELLS, &gx, &gy);
      a[nn*NUM_DIMS + 0] += gx/rho_p;
      a[nn*NUM_DIMS + 1] += gy/rho_p;

      // du/dt = -(P/rho) div(v), already divided by density on the mesh
      du[nn - species[SPECIES_GAS].offset] = cic_interp(du_mesh, ii, jj, dr, NUM_CELLS);
    }
    else {
      du[nn - species[SPECIES_GAS].offset] = 0.0;
    }
  }

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

#if !SHOCK_TUBE
  // initial condition variables
  double r2[NUM_DIMS];
  double dd;
  double vv;
#endif
  
  double *rho, *sub_rho, *rho_species;
  double *u, *du;
  double *press, *du_mesh, *gas_vel;
  double cs_max, v_max, e_thermal;
  double acc;
  int s;
  double *phi;
  double *green;
  double *gr_hat;

  double k_x, k_y;
  
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
  if(!press || !du_mesh || !gas_vel){
    printf("Out of memory: gas pressure fields not allocated.\n");
    exit(1);
  }
  // potential
  phi=(double *)malloc(NUM_CELLS*NUM_CELLS*sizeof(double));
  if(!phi){
    printf("Out of memory: phi not allocated.\n");
    exit(1);
  }
  // green's function
  green=(double *)malloc((NUM_CELLS)*(NUM_CELLS)*sizeof(double));
  if(!green){
    printf("Out of memory: green not allocated.\n");
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
  
  // compute gravitational green's function in k-space
  for(jj=0; jj<NUM_CELLS/2; jj++) {
    for(ii=0; ii<NUM_CELLS/2; ii++) {
      k_x = 2.0*M_PI*(ii+1)/(double)(NUM_CELLS);
      k_y = 2.0*M_PI*(jj+1)/(double)(NUM_CELLS);

      green[jj*NUM_CELLS + ii] = -(1.0/4.0)*1.0/(sin(k_x/2.0)*sin(k_x/2.0)+sin(k_y/2.0)*sin(k_y/2.0));
      gr_hat[jj*NUM_CELLS + ii] = green[jj*NUM_CELLS + ii];
      gr_hat[jj*NUM_CELLS + (NUM_CELLS-ii-1)] = green[jj*NUM_CELLS + ii];
      gr_hat[(NUM_CELLS-jj-1)*NUM_CELLS + (NUM_CELLS-ii-1)] = green[jj*NUM_CELLS + ii];      
      gr_hat[(NUM_CELLS-jj-1)*NUM_CELLS + ii] = green[jj*NUM_CELLS + ii];      
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
    compute_gas_pressure(rho_species, press, du_mesh, gas_vel, &cs_max, &v_max);

    for(nn=0;nn<NUM_THREADS;nn++){
      thread_data_array[nn].r=r;
      thread_data_array[nn].a=a;
      thread_data_array[nn].du=du;
      thread_data_array[nn].grid_fields=rho_species;
      thread_data_array[nn].press=press;
      thread_data_array[nn].du_mesh=du_mesh;
    }
    run_threads(hydro_interpolation_thread);

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

