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

// test problems, all 0: protostellar gas disk around a central star
// SHOCK_TUBE 1: Sod shock tube along x, gas only, no gravity, writes shocktube.dat
// GAS_CLOUD 1: polytropic self-gravitating gas cylinder in hydrostatic
// equilibrium, writes radial profiles to gascloud.dat
// DUSTY_BOX 1: uniform gas and dust streaming through each other, the drag
// test of Laibe and Price (2011), writes dustybox.dat. With DB_WAVE > 0 a
// sound wave in the dusty gas instead (DUSTYWAVE)
#ifndef SHOCK_TUBE
#define SHOCK_TUBE 0
#endif
#ifndef GAS_CLOUD
#define GAS_CLOUD 0
#endif
#ifndef DUSTY_BOX
#define DUSTY_BOX 0
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

// pure APIC transfers
#ifndef PIC_BLEND
#define PIC_BLEND 1.0
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

// FLIP blending without APIC: with APIC transfers the hot self-gravitating
// cloud slowly builds grid scale velocity noise, most of it in the affine
// momentum, that FLIP blending keeps small. The growth does not go away at
// smaller time steps, it is a spatial instability not yet understood
#ifndef APIC
#define APIC 0
#endif
#elif DUSTY_BOX
// uniform gas at rest and uniform dust moving through it with DB_V0 along x,
// no gravity and no pressure gradient, so only the drag acts: the center of
// mass velocity stays constant and the relative velocity decays as
// exp(-(1 + eps) t/t_s) with the dust to gas ratio eps, the kinetic energy
// lost heats the gas. Regular lattices of DB_PPC_X by DB_PPC_Y gas and
// dust particles per cell
//
// DB_WAVE > 0: DUSTYWAVE, gas and dust at rest on average carry a sound wave
// of DB_WAVE wavelengths along x, both displaced by the same amplitude
// DB_AMP and moving with the gas velocity of a wave travelling to +x. The
// linear solution of the coupled waves is computed by tests/dustywave.py
// from the Fourier amplitudes written to dustywave.dat
#define DB_PPC_X 2
#define DB_PPC_Y 2
#ifndef DB_WAVE
#define DB_WAVE 0
#endif
#if DB_WAVE
#define DB_CS 10.0
#define DB_V0 0.0
#define DB_AMP 1.0e-3
#define DB_K (2.0*M_PI*DB_WAVE/(double)(NUM_CELLS))
#else
#define DB_CS 1.0
#define DB_V0 1.0
#endif
#ifndef DB_EPS
#define DB_EPS 1.0
#endif
// constant dust stopping time
#ifndef DUST_TS
#define DUST_TS 0.05
#endif

#define NUM_DM_PARTICLES 0
#define NUM_GAS_PARTICLES (NUM_CELLS*NUM_CELLS*DB_PPC_X*DB_PPC_Y)
#define NUM_DUST_PARTICLES (NUM_CELLS*NUM_CELLS*DB_PPC_X*DB_PPC_Y)
#define GAS_MASS_FRACTION 1.0
#define DUST_TO_GAS DB_EPS
#define GRAVITY 0
#ifndef T_END
#if DB_WAVE
// one sound crossing of a wavelength
#define T_END (NUM_CELLS/(DB_WAVE*DB_CS))
#else
#define T_END 0.5
#endif
#endif

// pure APIC transfers
#ifndef PIC_BLEND
#define PIC_BLEND 1.0
#endif
#else
// star system formation: a self-gravitating gas disk around a protostar, no
// dark matter. The collisionless species holds a single particle, the star,
// which moves freely under the mesh gravity of the disk. The disk surface
// density is sigma ~ 1/r between SD_R_IN and SD_R_OUT with tanh edges of
// width SD_EDGE_IN and SD_EDGE_OUT. The wide outer taper keeps the pressure
// support of the outer edge small, a sharp hot edge loses it when it cools,
// falls inwards and piles up into a ring that fragments before the disk.
//
// Gravity is that of a thin disk in 3D (THIN_DISK_GRAVITY): the potential of
// Plummer softened point masses -G m/sqrt(r^2 + h^2), whose transform is
// phi_k = -2 pi G sigma_k exp(-k h)/k, so the star's pull falls off as 1/r^2
// and the rotation is Keplerian. It is solved on a zero padded mesh without
// periodic images, whose slowly decaying 1/r tides would impose the square
// symmetry of the box on the disk. With the 2D log potential of the other
// test problems a rotating disk cannot fragment, kappa^2 > 4 pi G sigma always.
//
// The gas is an ideal gas cooled on beta cooling times t_cool = SD_BETA/Omega
// towards the temperature that stellar irradiation keeps it at, a function of
// radius only: cs_floor ~ r^(-1/4) (T ~ r^(-1/2)), normalized to Toomre
// Q = cs kappa/(pi G sigma) = SD_Q_FLOOR in the middle of the disk. With
// sigma ~ 1/r the floor Q falls as r^(-3/4), so the outer disk can fragment
// while the inner disk stays marginal and the thin outer taper stays warm and
// stable. The gas starts at Q = SD_Q_INIT but not below the floor
// temperature. Fast cooling (SD_BETA of a few) lets Q fall below ~1,
// the disk breaks up into spiral arms and bound gas clumps collapse out of
// them; slow cooling heats the disk by gravitational turbulence instead
// (Gammie 2001). The rotation is put into equilibrium from the gravity and
// pressure forces the code computes for the initial positions (see
// set_equilibrium_rotation).
#define STAR_DISK 1
#ifndef SD_DISK_TO_STAR
#define SD_DISK_TO_STAR 0.2
#endif
#define SD_R_IN 16.0
#define SD_R_OUT 80.0
#define SD_EDGE_IN 3.0
#define SD_EDGE_OUT 12.0
#define SD_OUT_DT 4.0
#define SD_STAR_MASS (1.0/(1.0 + SD_DISK_TO_STAR))

// thin disk gravity with G = SD_G_FACTOR G_EFF and softening SD_SOFT cells,
// the large G keeps the inner disk orbit near the maximum time step
#define THIN_DISK_GRAVITY 1
#define SD_G_FACTOR 20.0
#define G_DISK (SD_G_FACTOR*G_EFF)
#define SD_SOFT 2.0
// circular speed of the star alone at the outer disk radius
#define SD_V_K sqrt(G_DISK*SD_STAR_MASS/SD_R_OUT)

// ideal gas with cooling
#define ISOTHERMAL 0
#define GAMMA (5.0/3.0)
#ifndef SD_Q_INIT
#define SD_Q_INIT 1.5
#endif
#ifndef SD_COOLING
#define SD_COOLING 1
#endif
#ifndef SD_BETA
#define SD_BETA 2.0
#endif
#ifndef SD_Q_FLOOR
#define SD_Q_FLOOR 0.5
#endif

// gas within SD_R_ACC cells of the star is accreted onto it
#define SD_R_ACC 8.0

// pure APIC transfers: FLIP blending toward the CIC mesh velocity exerts a
// net torque on the disk proportional to the blend fraction, APIC keeps the
// rotation in the affine velocity
#ifndef PIC_BLEND
#define PIC_BLEND 1.0
#endif

#define NUM_DM_PARTICLES 1
#define NUM_GAS_PARTICLES 1000000
#define GAS_MASS_FRACTION (SD_DISK_TO_STAR/(1.0 + SD_DISK_TO_STAR))

// dust in the disk (SD_DUST 1): superparticles placed like the gas and
// rotating with it, DUST_TO_GAS of the gas mass. Small grains in the Epstein
// regime of a thin disk have the Stokes number St = Omega t_s =
// pi a rho_s/(2 sigma_gas), independent of the temperature; DUST_ST is St of
// the grain size in the initial gas surface density in the middle of the disk
#ifndef SD_DUST
#define SD_DUST 0
#endif
#if SD_DUST
#define NUM_DUST_PARTICLES 500000
#endif
#ifndef DUST_ST
#define DUST_ST 1.0
#endif

// self-gravity on, run until the window is closed unless T_END is given
#define GRAVITY 1
#ifndef T_END
#define T_END 0.0
#endif
#endif

#ifndef STAR_DISK
#define STAR_DISK 0
#endif
#ifndef THIN_DISK_GRAVITY
#define THIN_DISK_GRAVITY 0
#endif

// dust superparticles, DUST_TO_GAS of the gas mass in total. A positive
// DUST_TS is a constant stopping time, otherwise the test problem sets it
#ifndef NUM_DUST_PARTICLES
#define NUM_DUST_PARTICLES 0
#endif
#ifndef DUST_TO_GAS
#define DUST_TO_GAS 0.01
#endif
#ifndef DUST_TS
#define DUST_TS 0.0
#endif

#define NUM_PARTICLES (NUM_DM_PARTICLES + NUM_GAS_PARTICLES + NUM_DUST_PARTICLES)

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
#define DUST_MASS ((NUM_DUST_PARTICLES) > 0 ? DUST_TO_GAS*GAS_MASS_FRACTION/(double)(NUM_DUST_PARTICLES) : 0.0)

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
// applied in converging flow only. It acts through the cell faces, measuring
// the compression with the velocity divergence at the face: the normal
// velocity jump across the face, a compact difference that also damps the
// odd-even velocity mode the central difference pressure force cannot see
// and APIC transfers do not filter, plus the tangential derivative of the
// tangential velocity, central differences averaged over the two cells. The
// normal jump alone sees compression in pure shear, in a Keplerian disk
// dvx/dx = -dvy/dy = 3/4 Omega sin(2 phi) per cell, and heats the disk on a
// fraction of an orbit where its scale height is unresolved. AV_FACE_DIV 0
// returns to the normal jump alone
#ifndef AV_FACE_DIV
#define AV_FACE_DIV 1
#endif
#ifndef Q_LIN
#define Q_LIN 0.5
#endif
#ifndef Q_QUAD
#define Q_QUAD 2.0
#endif

// FLIP/PIC blending: each step gas particle velocities relax by this fraction
// towards the mass weighted mesh velocity, 0 is pure FLIP (particles keep
// their own velocities), 1 is pure PIC (pure APIC with APIC transfers);
// damps sub-cell particle noise
#ifndef PIC_BLEND
#define PIC_BLEND 0.15
#endif

// APIC affine particle-in-cell transfers (Jiang et al. 2015): gas particles
// carry an affine velocity matrix B, the momentum deposit adds the particle's
// linear velocity field, m B grad(w), and the blend relaxes B towards the
// mesh value sum_c w V_c (x_c - x_p)^T by the same fraction as the velocity.
// Linear flows such as shear and rotation are then transferred without loss
// and the angular momentum sum m x_p x v_p + m (B_yx - B_xy) is conserved
// exactly for any blend fraction. Between blends B is not evolved, so it
// lags the flow for small blend fractions: use APIC with PIC_BLEND 1 (pure
// APIC) or close to it
#ifndef APIC
#define APIC 1
#endif

// the affine deposit m B grad(w) does not vanish with the CIC weight w, so a
// node that receives only a small weight from a lone particle gets a velocity
// mom/rho far beyond the flow. Where the sparsest node of a particle's stencil
// holds less than APIC_SPARSE particle masses its B is scaled down by
// rho_min/(APIC_SPARSE m). The B grad(w) terms of a particle sum to zero, so
// momentum and energy stay exact, only the spin of near-empty regions is lost
#ifndef APIC_SPARSE
#define APIC_SPARSE 4.0
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

// species table, dark matter first, gas after it, dust last
struct species_info species[NUM_SPECIES] = {
  [SPECIES_DM]   = {0, NUM_DM_PARTICLES, DM_MASS},
  [SPECIES_GAS]  = {NUM_DM_PARTICLES, NUM_GAS_PARTICLES, GAS_MASS},
  [SPECIES_DUST] = {NUM_DM_PARTICLES + NUM_GAS_PARTICLES, NUM_DUST_PARTICLES, DUST_MASS},
};

// mass of particle nn
static inline double particle_mass(int nn){
  if(nn < species[SPECIES_GAS].offset) {
    return species[SPECIES_DM].mass;
  }
  if(nn < species[SPECIES_DUST].offset) {
    return species[SPECIES_GAS].mass;
  }
  return species[SPECIES_DUST].mass;
}

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
  double *r, *v, *u, *aff, *sub_rho, *rho, *rho_dm;
  double *e_gas, *mx_gas, *my_gas, *ax_gas, *ay_gas, *mx_dust, *my_dust, *ax_dust, *ay_dust;
  const unsigned char *sink, *dsink;
  double dr[2];
  double w[4], gw[4][2], d[4][2];
  double mass;
  const double *bp;
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
  aff=my_data->aff;
  sink=my_data->sink;
  dsink=my_data->dsink;
  sub_rho=my_data->rho;

  thread_id=my_data->thread_id;

  // gas energy and momentum subfields of this thread
  e_gas=&sub_rho[(FIELD_GAS_ENERGY*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];
  mx_gas=&sub_rho[(FIELD_GAS_MOM_X*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];
  my_gas=&sub_rho[(FIELD_GAS_MOM_Y*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];
  ax_gas=&sub_rho[(FIELD_GAS_AFF_X*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];
  ay_gas=&sub_rho[(FIELD_GAS_AFF_Y*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];
  mx_dust=&sub_rho[(FIELD_DUST_MOM_X*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];
  my_dust=&sub_rho[(FIELD_DUST_MOM_Y*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];
  ax_dust=&sub_rho[(FIELD_DUST_AFF_X*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];
  ay_dust=&sub_rho[(FIELD_DUST_AFF_Y*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];
  rho_dm=&sub_rho[(SPECIES_DM*NUM_THREADS + thread_id)*NUM_CELLS*NUM_CELLS];

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

      cic_stencil(ii, jj, dr, NUM_CELLS, c, w, gw, d);
      for(kk=0; kk<4; kk++) {
	w[kk] *= mass;
      }

      // gas and dust accreted by a sink are collisionless mass of the sink
      if((s == SPECIES_GAS && sink[nn - species[SPECIES_GAS].offset]) ||
	 (s == SPECIES_DUST && dsink[nn - species[SPECIES_DUST].offset])) {
	for(kk=0; kk<4; kk++) {
	  rho_dm[c[kk]] += w[kk];
	}
	continue;
      }

      for(kk=0; kk<4; kk++) {
	rho[c[kk]] += w[kk];
      }

      // gas also deposits internal energy and momentum, with APIC the
      // momentum of the particle's affine velocity field m B grad(w) goes
      // to separate fields
      if(s == SPECIES_GAS) {
	bp = &aff[4*(nn - species[SPECIES_GAS].offset)];
	for(kk=0; kk<4; kk++) {
	  e_gas[c[kk]] += w[kk]*u[nn - species[SPECIES_GAS].offset];
	  mx_gas[c[kk]] += w[kk]*v[nn*NUM_DIMS + 0];
	  my_gas[c[kk]] += w[kk]*v[nn*NUM_DIMS + 1];
	  if(APIC) {
	    ax_gas[c[kk]] += mass*(bp[0]*gw[kk][0] + bp[1]*gw[kk][1]);
	    ay_gas[c[kk]] += mass*(bp[2]*gw[kk][0] + bp[3]*gw[kk][1]);
	  }
	}
      }

      // dust momentum for the drag, with its APIC affine momentum
      if(s == SPECIES_DUST) {
	bp = &my_data->daff[4*(nn - species[SPECIES_DUST].offset)];
	for(kk=0; kk<4; kk++) {
	  mx_dust[c[kk]] += w[kk]*v[nn*NUM_DIMS + 0];
	  my_dust[c[kk]] += w[kk]*v[nn*NUM_DIMS + 1];
	  ax_dust[c[kk]] += mass*(bp[0]*gw[kk][0] + bp[1]*gw[kk][1]);
	  ay_dust[c[kk]] += mass*(bp[2]*gw[kk][0] + bp[3]*gw[kk][1]);
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

// artificial viscous pressure of the face between cells i and j with
// densities rho_i, rho_j, thermal pressures p_i, p_j and face velocity
// divergence dv (times the cell size), nonzero in compression only,
// symmetric under exchanging the cells so both cells see the same value
static inline double face_q(double rho_i, double rho_j, double p_i, double p_j, double dv){
  double cs;

  if(dv >= 0.0) {
    return 0.0;
  }
#if ISOTHERMAL
  cs = CS_ISO;
  (void)p_i;
  (void)p_j;
#else
  cs = sqrt(fmax(GAMMA*(p_i + p_j)/(rho_i + rho_j), 0.0));
#endif
  return 0.5*(rho_i + rho_j)*(Q_LIN*cs*(-dv) + Q_QUAD*dv*dv);
}

// compute gas pressure including artificial viscosity, the face pressure
// acceleration and the matching specific internal energy rate on the mesh,
// returns maximum sound speed and velocity and the residual of the discrete
// work balance, which is zero up to roundoff
void compute_gas_pressure(const double *grid_fields, double *press, double *du_mesh, double *vel, double *acc, double *cs_max, double *v_max, double *p_res){
  int ii, jj, nn;
  int nxp, nxm, nyp, nym;
  double rho, cs2, vabs;
  double l_xm, l_xp, l_ym, l_yp;
  double q_xm, q_xp, q_ym, q_yp;
  double dv_xm, dv_xp, dv_ym, dv_yp;
  // central differences dvx/dx and dvy/dy per cell for the tangential part
  // of the face divergence
  static double dxvx[NUM_CELLS*NUM_CELLS], dyvy[NUM_CELLS*NUM_CELLS];
  double tg_xm, tg_xp, tg_ym, tg_yp;

  const double *rho_gas = &grid_fields[FIELD_RHO_GAS*NUM_CELLS*NUM_CELLS];
  const double *e_gas = &grid_fields[FIELD_GAS_ENERGY*NUM_CELLS*NUM_CELLS];
  const double *mx_gas = &grid_fields[FIELD_GAS_MOM_X*NUM_CELLS*NUM_CELLS];
  const double *my_gas = &grid_fields[FIELD_GAS_MOM_Y*NUM_CELLS*NUM_CELLS];
  const double *affx_gas = &grid_fields[FIELD_GAS_AFF_X*NUM_CELLS*NUM_CELLS];
  const double *affy_gas = &grid_fields[FIELD_GAS_AFF_Y*NUM_CELLS*NUM_CELLS];
  // full mesh velocity (mom + affine)/rho for the blend and the CFL limit,
  // and particle velocity part mom/rho for the viscosity and the face work
  double *vfx = &vel[0];
  double *vfy = &vel[NUM_CELLS*NUM_CELLS];
  double *vx = &vel[2*NUM_CELLS*NUM_CELLS];
  double *vy = &vel[3*NUM_CELLS*NUM_CELLS];
  double *ax = &acc[0];
  double *ay = &acc[NUM_CELLS*NUM_CELLS];

  *cs_max = 0.0;
  *v_max = 0.0;

  // mass weighted mesh velocity and thermal pressure. The kick changes the
  // particle velocities only, so the work of the face forces is taken with
  // the particle velocity part of the mesh velocity, the APIC affine
  // momentum receives no force and does no work
  for(nn=0; nn<NUM_CELLS*NUM_CELLS; nn++) {
    rho = rho_gas[nn];
    if(rho > 0.0) {
      vx[nn] = mx_gas[nn]/rho;
      vy[nn] = my_gas[nn]/rho;
      vfx[nn] = (mx_gas[nn] + affx_gas[nn])/rho;
      vfy[nn] = (my_gas[nn] + affy_gas[nn])/rho;
    }
    else {
      vx[nn] = 0.0;
      vy[nn] = 0.0;
      vfx[nn] = 0.0;
      vfy[nn] = 0.0;
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

  // maximum sound speed and velocity of the gas
  for(nn=0; nn<NUM_CELLS*NUM_CELLS; nn++) {
    rho = rho_gas[nn];

    if(rho <= RHO_GAS_FLOOR) {
      continue;
    }

#if ISOTHERMAL
    cs2 = CS_ISO*CS_ISO;
#else
    cs2 = GAMMA*press[nn]/rho;
#endif

    if(cs2 > *cs_max) {
      *cs_max = cs2;
    }
    vabs = vfx[nn]*vfx[nn] + vfy[nn]*vfy[nn];
    if(vabs > *v_max) {
      *v_max = vabs;
    }
  }

  // tangential derivatives of the particle velocity part of the mesh velocity
  for(jj=0; jj<NUM_CELLS; jj++) {
    for(ii=0; ii<NUM_CELLS; ii++) {
      nn = jj*NUM_CELLS + ii;
#if AV_FACE_DIV
      dxvx[nn] = 0.5*(vx[jj*NUM_CELLS + modulo(ii+1, NUM_CELLS)] - vx[jj*NUM_CELLS + modulo(ii-1, NUM_CELLS)]);
      dyvy[nn] = 0.5*(vy[modulo(jj+1, NUM_CELLS)*NUM_CELLS + ii] - vy[modulo(jj-1, NUM_CELLS)*NUM_CELLS + ii]);
#else
      dxvx[nn] = 0.0;
      dyvy[nn] = 0.0;
#endif
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

      dv_xm = vx[nn] - vx[nxm];
      dv_xp = vx[nxp] - vx[nn];
      dv_ym = vy[nn] - vy[nym];
      dv_yp = vy[nyp] - vy[nn];

      // face divergence: normal jump plus the tangential derivative
      // averaged over the two cells of the face
      tg_xm = 0.5*(dyvy[nxm] + dyvy[nn]);
      tg_xp = 0.5*(dyvy[nn] + dyvy[nxp]);
      tg_ym = 0.5*(dxvx[nym] + dxvx[nn]);
      tg_yp = 0.5*(dxvx[nn] + dxvx[nyp]);

      q_xm = face_q(rho_gas[nxm], rho, press[nxm], press[nn], dv_xm + tg_xm);
      q_xp = face_q(rho, rho_gas[nxp], press[nn], press[nxp], dv_xp + tg_xp);
      q_ym = face_q(rho_gas[nym], rho, press[nym], press[nn], dv_ym + tg_ym);
      q_yp = face_q(rho, rho_gas[nyp], press[nn], press[nyp], dv_yp + tg_yp);

      ax[nn] = (l_xm*(0.5*(press[nxm] + press[nn]) + q_xm) - l_xp*(0.5*(press[nn] + press[nxp]) + q_xp))/rho;
      ay[nn] = (l_ym*(0.5*(press[nym] + press[nn]) + q_ym) - l_yp*(0.5*(press[nn] + press[nyp]) + q_yp))/rho;

      // thermal face work split in proportion to the cell pressures, the
      // viscous face work equally between the two cells
      du_mesh[nn] = -0.5*(l_xp*(press[nn] + q_xp)*dv_xp + l_xm*(press[nn] + q_xm)*dv_xm +
			  l_yp*(press[nn] + q_yp)*dv_yp + l_ym*(press[nn] + q_ym)*dv_ym)/rho;
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
  double vpx, vpy, wx, wy, var, v2, h, bc;
  double gw[4][2], d[4][2], bg[4], cg[4];
  double *aff, *bp;
  const double *rho_gas;
  double rho_min, sparse;
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
  aff=my_data->aff;
  rho_gas=&my_data->grid_fields[FIELD_RHO_GAS*NUM_CELLS*NUM_CELLS];

  thread_id=my_data->thread_id;

  thread_slice(species[SPECIES_GAS].offset, species[SPECIES_GAS].count, thread_id, NUM_THREADS, &lo, &hi);

  for(nn=lo; nn<hi; nn++) {
    ii = floor(r[nn*NUM_DIMS + 0]);
    jj = floor(r[nn*NUM_DIMS + 1]);

    dr[0] = r[nn*NUM_DIMS + 0] - (double)(ii);
    dr[1] = r[nn*NUM_DIMS + 1] - (double)(jj);

    gg = nn - species[SPECIES_GAS].offset;

    // accreted gas only feels gravity
    if(my_data->sink[gg]) {
      du[gg] = 0.0;
      continue;
    }

    // a = f/rho from the face pressure forces, divided by density on the mesh
    a[nn*NUM_DIMS + 0] += cic_interp(ax, ii, jj, dr, NUM_CELLS);
    a[nn*NUM_DIMS + 1] += cic_interp(ay, ii, jj, dr, NUM_CELLS);

    // du/dt from the face work, divided by density on the mesh
    du[gg] = cic_interp(du_mesh, ii, jj, dr, NUM_CELLS);

    if(PIC_BLEND > 0.0) {
      cic_stencil(ii, jj, dr, NUM_CELLS, c, w, gw, d);
      bp = &aff[4*gg];

      // mass weighted mesh velocity V at the particle, its CIC weighted
      // variance over the four cells, the mesh velocity gradient
      // C_ab = sum_c V_c,a dw_c/dx_b and the APIC matrix of the mesh
      // velocity B^g_ab = sum_c w_c V_c,a (x_c - x_p)_b
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
      bc = 0.0;
      if(APIC) {
	bg[0] = bg[1] = bg[2] = bg[3] = 0.0;
	cg[0] = cg[1] = cg[2] = cg[3] = 0.0;
	for(kk=0; kk<4; kk++) {
	  bg[0] += w[kk]*vx[c[kk]]*d[kk][0];
	  bg[1] += w[kk]*vx[c[kk]]*d[kk][1];
	  bg[2] += w[kk]*vy[c[kk]]*d[kk][0];
	  bg[3] += w[kk]*vy[c[kk]]*d[kk][1];
	  cg[0] += vx[c[kk]]*gw[kk][0];
	  cg[1] += vx[c[kk]]*gw[kk][1];
	  cg[2] += vy[c[kk]]*gw[kk][0];
	  cg[3] += vy[c[kk]]*gw[kk][1];
	}
	// B : C with the B that was deposited
	bc = bp[0]*cg[0] + bp[1]*cg[1] + bp[2]*cg[2] + bp[3]*cg[3];
	rho_min = rho_gas[c[0]];
	for(kk=1; kk<4; kk++) {
	  rho_min = fmin(rho_min, rho_gas[c[kk]]);
	}
	sparse = fmin(1.0, rho_min/(APIC_SPARSE*species[SPECIES_GAS].mass));
	for(kk=0; kk<4; kk++) {
	  bp[kk] += PIC_BLEND*(bg[kk] - bp[kk]);
	  bp[kk] *= sparse;
	}
      }

      // velocity relative to the mesh velocity, blended by a = PIC_BLEND
      wx = v[nn*NUM_DIMS + 0] - vpx;
      wy = v[nn*NUM_DIMS + 1] - vpy;

      v2 = v[nn*NUM_DIMS + 0]*v[nn*NUM_DIMS + 0] + v[nn*NUM_DIMS + 1]*v[nn*NUM_DIMS + 1];
      v[nn*NUM_DIMS + 0] -= PIC_BLEND*wx;
      v[nn*NUM_DIMS + 1] -= PIC_BLEND*wy;
      flip_loss += 0.5*(v2 - v[nn*NUM_DIMS + 0]*v[nn*NUM_DIMS + 0] - v[nn*NUM_DIMS + 1]*v[nn*NUM_DIMS + 1]);

      // the removed kinetic energy becomes heat. Per particle the blend removes
      // a(2-a)/2 |w|^2 + a V.w. Because V_c = mom_c/rho_c and the deposit
      // carried m (w v_p + B grad(w)), summed over particles the second term
      // is a sum_p m (var_p - B:C_p), so heating each particle by
      // a(2-a)/2 |w|^2 + a (var_p - B:C_p) returns exactly the kinetic energy
      // removed with a local and Galilean invariant heat. Without APIC B = 0
      // and the heat is non-negative, with APIC var_p - B:C_p vanishes for a
      // linear velocity field, so resolved shear is not dissipated
      if(PIC_BLEND_HEAT && !ISOTHERMAL) {
	h = 0.5*PIC_BLEND*(2.0 - PIC_BLEND)*(wx*wx + wy*wy) + PIC_BLEND*(var - bc);
	u[gg] = fmax(u[gg] + h, U_MIN);
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

#if STAR_DISK
// star disk surface density profile without normalization
static double star_disk_sigma(double rr){
  return 0.25*(1.0 + tanh((rr - SD_R_IN)/SD_EDGE_IN))*(1.0 - tanh((rr - SD_R_OUT)/SD_EDGE_OUT))/rr;
}

#define SD_R_MAX (SD_R_OUT + 6.0*SD_EDGE_OUT)
#define SD_NUM_TAB 8192

// initial and irradiation floor specific internal energy and orbital
// frequency on the radius grid
// rr = SD_R_MAX k/SD_NUM_TAB, k = 0..SD_NUM_TAB, set by init_star_disk
static double sd_u_tab[SD_NUM_TAB + 1];
static double sd_ufloor_tab[SD_NUM_TAB + 1];
static double sd_omega_tab[SD_NUM_TAB + 1];

// gas surface density of the initial profile in the middle of the disk, the
// reference of DUST_ST, set by init_star_disk
static double sd_sigma_ref;

// linear interpolation of a radius grid table, constant beyond its ends
static double sd_lookup(const double *tab, double rr){
  double f = rr/SD_R_MAX*(double)(SD_NUM_TAB);
  int kk = (int)f;

  if(kk >= SD_NUM_TAB) {
    return tab[SD_NUM_TAB];
  }
  f -= (double)(kk);
  return (1.0 - f)*tab[kk] + f*tab[kk+1];
}

// sound speed of the initial temperature profile at radius rr
static double sd_cs_init(double rr){
  return sqrt(GAMMA*(GAMMA - 1.0)*sd_lookup(sd_u_tab, rr));
}

// uniform random number in [0, 1) from a 64 bit xorshift* generator
static double sd_uniform(unsigned long long *state){
  *state ^= *state >> 12;
  *state ^= *state << 25;
  *state ^= *state >> 27;
  return (double)((*state*2685821657736338717ULL) >> 11)*(1.0/9007199254740992.0);
}

// star at rest in the box center, gas at radii from the inverse of the
// tabulated enclosed mass of the surface density profile and random
// azimuths: the azimuthal noise seeds the gravitational instability at all
// angular scales, a smooth (golden angle) pattern would leave the weak four
// fold anisotropy of the Cartesian mesh as its only seed. Velocities are set
// later by set_equilibrium_rotation. The temperature is
// set for Toomre Q = SD_Q_INIT from the Keplerian estimates
// Omega^2 = G M(<r)/r^3 and kappa^2 = Omega^2 + 2 pi G sigma/r of the
// enclosed star and disk mass, softened like the mesh gravity, but not
// below the irradiation floor, cs_floor ~ r^(-1/4) with Q = SD_Q_FLOOR at
// the middle of the disk
void init_star_disk(double *r, double *v, double *a, double *u, double *du){
  int nn, gg, kk;
  double rr, q, th, f, sigma, menc, s3, om2, kap2, cs, cs_ref, r_ref;
  unsigned long long seed = 0x9e3779b97f4a7c15ULL;
  static double rad_tab[SD_NUM_TAB + 1], cdf_tab[SD_NUM_TAB + 1];
  // pi G sigma/kappa, the sound speed of Q = 1
  static double qcs_tab[SD_NUM_TAB + 1];

  // enclosed mass table, sigma ~ 1/r vanishes towards r = 0 through the taper
  rad_tab[0] = 0.0;
  cdf_tab[0] = 0.0;
  for(kk=1; kk<=SD_NUM_TAB; kk++) {
    rad_tab[kk] = SD_R_MAX*(double)(kk)/(double)(SD_NUM_TAB);
    rr = 0.5*(rad_tab[kk] + rad_tab[kk-1]);
    cdf_tab[kk] = cdf_tab[kk-1] + 2.0*M_PI*rr*star_disk_sigma(rr)*(rad_tab[kk] - rad_tab[kk-1]);
  }

  // temperature and orbital frequency tables
  for(kk=0; kk<=SD_NUM_TAB; kk++) {
    rr = rad_tab[kk];
    sigma = GAS_MASS_FRACTION*(rr > 0.0 ? star_disk_sigma(rr) : 0.0)/cdf_tab[SD_NUM_TAB];
    menc = SD_STAR_MASS + GAS_MASS_FRACTION*cdf_tab[kk]/cdf_tab[SD_NUM_TAB];
    s3 = pow(rr*rr + SD_SOFT*SD_SOFT, 1.5);
    om2 = G_DISK*menc/s3;
    kap2 = om2 + 2.0*M_PI*G_DISK*sigma*rr*rr/s3;
    qcs_tab[kk] = M_PI*G_DISK*sigma/sqrt(kap2);
    sd_omega_tab[kk] = sqrt(om2);
  }
  r_ref = 0.5*(SD_R_IN + SD_R_OUT);
  cs_ref = SD_Q_FLOOR*sd_lookup(qcs_tab, r_ref);
  sd_sigma_ref = GAS_MASS_FRACTION*star_disk_sigma(r_ref)/cdf_tab[SD_NUM_TAB];
  for(kk=0; kk<=SD_NUM_TAB; kk++) {
    // the floor is held constant inside the accretion radius
    cs = cs_ref*pow(fmax(rad_tab[kk], SD_R_ACC)/r_ref, -0.25);
    sd_ufloor_tab[kk] = cs*cs/(GAMMA*(GAMMA - 1.0));
    cs = SD_Q_INIT*qcs_tab[kk];
    sd_u_tab[kk] = fmax(fmax(cs*cs/(GAMMA*(GAMMA - 1.0)), sd_ufloor_tab[kk]), U_MIN);
  }

  for(nn=0; nn<NUM_PARTICLES; nn++) {
    v[nn*NUM_DIMS + 0] = 0.0;
    v[nn*NUM_DIMS + 1] = 0.0;
    a[nn*NUM_DIMS + 0] = 0.0;
    a[nn*NUM_DIMS + 1] = 0.0;
  }

  // the star
  for(nn=species[SPECIES_DM].offset; nn<species[SPECIES_DM].offset + species[SPECIES_DM].count; nn++) {
    r[nn*NUM_DIMS + 0] = 0.5*NUM_CELLS;
    r[nn*NUM_DIMS + 1] = 0.5*NUM_CELLS;
  }

  kk = 1;
  for(gg=0; gg<NUM_GAS_PARTICLES; gg++) {
    nn = species[SPECIES_GAS].offset + gg;
    q = ((double)(gg) + 0.5)/(double)(NUM_GAS_PARTICLES)*cdf_tab[SD_NUM_TAB];
    while(cdf_tab[kk] < q) {
      kk++;
    }
    f = (q - cdf_tab[kk-1])/(cdf_tab[kk] - cdf_tab[kk-1]);
    rr = rad_tab[kk-1] + f*(rad_tab[kk] - rad_tab[kk-1]);
    th = 2.0*M_PI*sd_uniform(&seed);

    r[nn*NUM_DIMS + 0] = 0.5*NUM_CELLS + rr*cos(th);
    r[nn*NUM_DIMS + 1] = 0.5*NUM_CELLS + rr*sin(th);
    u[gg] = sd_lookup(sd_u_tab, rr);
    du[gg] = 0.0;
  }

  // dust follows the gas surface density with its own random azimuths
  kk = 1;
  seed = 0x2545f4914f6cdd1dULL;
  for(gg=0; gg<NUM_DUST_PARTICLES; gg++) {
    nn = species[SPECIES_DUST].offset + gg;
    q = ((double)(gg) + 0.5)/(double)(NUM_DUST_PARTICLES)*cdf_tab[SD_NUM_TAB];
    while(cdf_tab[kk] < q) {
      kk++;
    }
    f = (q - cdf_tab[kk-1])/(cdf_tab[kk] - cdf_tab[kk-1]);
    rr = rad_tab[kk-1] + f*(rad_tab[kk] - rad_tab[kk-1]);
    th = 2.0*M_PI*sd_uniform(&seed);

    r[nn*NUM_DIMS + 0] = 0.5*NUM_CELLS + rr*cos(th);
    r[nn*NUM_DIMS + 1] = 0.5*NUM_CELLS + rr*sin(th);
  }
}

// beta cooling of the gas towards the floor temperature over the cooling
// time t_cool = SD_BETA/Omega at the particle's distance from the star,
// u_floor + (u - u_floor) exp(-dt/t_cool) solved exactly, returns the
// radiated energy (negative where the floor heats the gas)
double cool_star_disk(const double *r, double *u, const unsigned char *sink, double dt){
  int nn, gg;
  int ns = species[SPECIES_DM].offset;
  double cx = r[ns*NUM_DIMS + 0], cy = r[ns*NUM_DIMS + 1];
  double dx, dy, rr, u_floor, u_new, e_rad = 0.0;

  for(gg=0; gg<NUM_GAS_PARTICLES; gg++) {
    nn = species[SPECIES_GAS].offset + gg;
    if(sink[gg]) {
      continue;
    }
    dx = r[nn*NUM_DIMS + 0] - cx;
    dy = r[nn*NUM_DIMS + 1] - cy;
    rr = sqrt(dx*dx + dy*dy);
    u_floor = sd_lookup(sd_ufloor_tab, rr);
    u_new = u_floor + (u[gg] - u_floor)*exp(-dt*sd_lookup(sd_omega_tab, rr)/SD_BETA);
    u_new = fmax(u_new, U_MIN);
    e_rad += u[gg] - u_new;
    u[gg] = u_new;
  }
  return species[SPECIES_GAS].mass*e_rad;
}

// accrete gas and dust within SD_R_ACC of the star: the particle is moved
// onto the star and the star and all particles it accreted take the common
// velocity that conserves their momentum. Accreted particles then sit exactly
// on the star, receive the same gravitational acceleration and stay with it.
// Returns the number of particles accreted so far
int accrete_star_disk(double *r, double *v, double *aff, unsigned char *sink, unsigned char *dsink, double *daff){
  int nn, gg, kk, s;
  int ns = species[SPECIES_DM].offset;
  double cx = r[ns*NUM_DIMS + 0], cy = r[ns*NUM_DIMS + 1];
  double dx, dy, m;
  double px, py, mg;
  int n_acc = 0, n_new = 0;
  unsigned char *flag;

  mg = species[SPECIES_DM].mass;
  for(s=SPECIES_GAS; s<=SPECIES_DUST; s++) {
    flag = s == SPECIES_GAS ? sink : dsink;
    for(gg=0; gg<species[s].count; gg++) {
      n_acc += flag[gg];
      mg += flag[gg]*species[s].mass;
    }
  }
  px = mg*v[ns*NUM_DIMS + 0];
  py = mg*v[ns*NUM_DIMS + 1];

  for(s=SPECIES_GAS; s<=SPECIES_DUST; s++) {
    flag = s == SPECIES_GAS ? sink : dsink;
    m = species[s].mass;
    for(gg=0; gg<species[s].count; gg++) {
      if(flag[gg]) {
	continue;
      }
      nn = species[s].offset + gg;
      dx = r[nn*NUM_DIMS + 0] - cx;
      dy = r[nn*NUM_DIMS + 1] - cy;
      // the box is periodic but the star stays far from the boundary
      if(dx*dx + dy*dy < SD_R_ACC*SD_R_ACC) {
	flag[gg] = 1;
	px += m*v[nn*NUM_DIMS + 0];
	py += m*v[nn*NUM_DIMS + 1];
	mg += m;
	for(kk=0; kk<4; kk++) {
	  (s == SPECIES_GAS ? aff : daff)[4*gg + kk] = 0.0;
	}
	n_new++;
      }
    }
  }

  if(n_new > 0) {
    v[ns*NUM_DIMS + 0] = px/mg;
    v[ns*NUM_DIMS + 1] = py/mg;
    for(s=SPECIES_GAS; s<=SPECIES_DUST; s++) {
      flag = s == SPECIES_GAS ? sink : dsink;
      for(gg=0; gg<species[s].count; gg++) {
	if(flag[gg]) {
	  nn = species[s].offset + gg;
	  r[nn*NUM_DIMS + 0] = cx;
	  r[nn*NUM_DIMS + 1] = cy;
	  v[nn*NUM_DIMS + 0] = v[ns*NUM_DIMS + 0];
	  v[nn*NUM_DIMS + 1] = v[ns*NUM_DIMS + 1];
	}
      }
    }
  }
  return n_acc + n_new;
}

// append the gas surface density map as NUM_CELLS^2 floats, one frame per
// profile block of stardisk.dat
void write_star_disk_map(FILE *fp, const double *rho_gas){
  static float frame[NUM_CELLS*NUM_CELLS];
  int nn;

  for(nn=0; nn<NUM_CELLS*NUM_CELLS; nn++) {
    frame[nn] = (float)rho_gas[nn];
  }
  fwrite(frame, sizeof(float), NUM_CELLS*NUM_CELLS, fp);
  fflush(fp);
}

// bins of half a cell out to SD_R_MAX, as an integer constant expression
#define SD_NUM_BINS (2*((int)SD_R_OUT + 6*(int)SD_EDGE_OUT) + 2)

// circular velocities in equilibrium with the accelerations a of the initial
// positions: the gas acceleration towards the star averaged over radial bins
// of half a cell gives v_phi^2 = -r <a_r>, interpolated linearly between bin
// centers. This includes the softened mesh gravity of the star, the disk self
// gravity, the periodic box and the gas pressure. Dust takes the rotation
// of the gas around it, so it feels no drag at first, then sheds the
// pressure support it lacks and drifts. Net momentum is removed so
// the center of mass stays at rest, returns the maximum speed
double set_equilibrium_rotation(double *r, double *v, const double *a){
  int nn, bb, b0;
  double dx, dy, rr, ar, vphi2, f, vmax = 0.0;
  double px = 0.0, py = 0.0, mtot = 0.0, m;
  double cx = r[species[SPECIES_DM].offset*NUM_DIMS + 0];
  double cy = r[species[SPECIES_DM].offset*NUM_DIMS + 1];
  static double sum[SD_NUM_BINS];
  static int cnt[SD_NUM_BINS];
  static double prof[SD_NUM_BINS];

  for(bb=0; bb<SD_NUM_BINS; bb++) {
    sum[bb] = 0.0;
    cnt[bb] = 0;
  }
  for(nn=species[SPECIES_GAS].offset; nn<species[SPECIES_GAS].offset + species[SPECIES_GAS].count; nn++) {
    dx = r[nn*NUM_DIMS + 0] - cx;
    dy = r[nn*NUM_DIMS + 1] - cy;
    rr = sqrt(dx*dx + dy*dy);
    bb = (int)(2.0*rr);
    if(rr > 0.0 && bb < SD_NUM_BINS) {
      ar = (a[nn*NUM_DIMS + 0]*dx + a[nn*NUM_DIMS + 1]*dy)/rr;
      sum[bb] += -ar*rr;
      cnt[bb]++;
    }
  }
  for(bb=0; bb<SD_NUM_BINS; bb++) {
    prof[bb] = cnt[bb] > 0 ? sum[bb]/(double)(cnt[bb]) : -1.0;
  }

  // dust, placed after the gas, takes the rotation of the gas around it
  for(nn=species[SPECIES_GAS].offset; nn<species[SPECIES_DUST].offset + species[SPECIES_DUST].count; nn++) {
    dx = r[nn*NUM_DIMS + 0] - cx;
    dy = r[nn*NUM_DIMS + 1] - cy;
    rr = sqrt(dx*dx + dy*dy);

    // bin centers are at (bb + 1/2)/2
    f = 2.0*rr - 0.5;
    b0 = (int)floor(f);
    f -= (double)(b0);
    if(b0 >= 0 && b0 + 1 < SD_NUM_BINS && prof[b0] >= 0.0 && prof[b0+1] >= 0.0) {
      vphi2 = (1.0 - f)*prof[b0] + f*prof[b0+1];
    }
    else {
      // too few particles around, use the particle's own acceleration
      vphi2 = rr > 0.0 ? -(a[nn*NUM_DIMS + 0]*dx + a[nn*NUM_DIMS + 1]*dy) : 0.0;
    }
    vphi2 = fmax(vphi2, 0.0);

    if(rr > 0.0) {
      v[nn*NUM_DIMS + 0] = -dy/rr*sqrt(vphi2);
      v[nn*NUM_DIMS + 1] = dx/rr*sqrt(vphi2);
    }
  }

  // remove the net momentum of star and gas together
  for(nn=0; nn<NUM_PARTICLES; nn++) {
    m = particle_mass(nn);
    px += m*v[nn*NUM_DIMS + 0];
    py += m*v[nn*NUM_DIMS + 1];
    mtot += m;
  }
  for(nn=0; nn<NUM_PARTICLES; nn++) {
    v[nn*NUM_DIMS + 0] -= px/mtot;
    v[nn*NUM_DIMS + 1] -= py/mtot;
    vmax = fmax(vmax, sqrt(v[nn*NUM_DIMS + 0]*v[nn*NUM_DIMS + 0] + v[nn*NUM_DIMS + 1]*v[nn*NUM_DIMS + 1]));
  }

  return vmax;
}

// radial profile of the gas disk around the star in bins of one cell:
// surface density, mean azimuthal and radial velocity relative to the star,
// the rms radial velocity and the mass weighted sound speed, header with the
// star state, the total angular momentum about the center of mass including
// the APIC spin, the gas thermal energy and the energy radiated so far. With
// dust three more columns: dust surface density, mean azimuthal and radial
// dust velocity relative to the star, and the accreted dust mass in the header
void write_star_disk_profile(FILE *fp, const double *r, const double *v, const double *aff, const double *u,
			     const unsigned char *sink, const unsigned char *dsink, double e_rad, double t){
  int nn, bb;
  int ns = species[SPECIES_DM].offset;
  double cx = r[ns*NUM_DIMS + 0], cy = r[ns*NUM_DIMS + 1];
  double svx = v[ns*NUM_DIMS + 0], svy = v[ns*NUM_DIMS + 1];
  double dx, dy, dvx, dvy, rr, m, mtot = 0.0;
  double gx = 0.0, gy = 0.0, gvx = 0.0, gvy = 0.0, lz = 0.0, e_th = 0.0, m_acc = 0.0, m_acc_dust = 0.0;
  enum { NB = (int)SD_R_OUT + 6*(int)SD_EDGE_OUT + 1 };
  double mass[NB], vphi[NB], vr[NB], vr2[NB], csb[NB];
  double dmass[NB], dvphi[NB], dvr[NB];

  for(bb=0; bb<NB; bb++) {
    mass[bb] = vphi[bb] = vr[bb] = vr2[bb] = csb[bb] = 0.0;
    dmass[bb] = dvphi[bb] = dvr[bb] = 0.0;
  }

  // center of mass position and velocity
  for(nn=0; nn<NUM_PARTICLES; nn++) {
    m = particle_mass(nn);
    gx += m*r[nn*NUM_DIMS + 0];
    gy += m*r[nn*NUM_DIMS + 1];
    gvx += m*v[nn*NUM_DIMS + 0];
    gvy += m*v[nn*NUM_DIMS + 1];
    mtot += m;
  }
  gx /= mtot; gy /= mtot; gvx /= mtot; gvy /= mtot;

  for(nn=0; nn<NUM_PARTICLES; nn++) {
    m = particle_mass(nn);
    lz += m*((r[nn*NUM_DIMS + 0] - gx)*(v[nn*NUM_DIMS + 1] - gvy) - (r[nn*NUM_DIMS + 1] - gy)*(v[nn*NUM_DIMS + 0] - gvx));
    if(nn >= species[SPECIES_DUST].offset) {
      if(dsink[nn - species[SPECIES_DUST].offset]) {
	m_acc_dust += m;
	continue;
      }
      dx = r[nn*NUM_DIMS + 0] - cx;
      dy = r[nn*NUM_DIMS + 1] - cy;
      dvx = v[nn*NUM_DIMS + 0] - svx;
      dvy = v[nn*NUM_DIMS + 1] - svy;
      rr = sqrt(dx*dx + dy*dy);
      bb = (int)(rr);
      if(rr > 0.0 && bb < NB) {
	dmass[bb] += m;
	dvphi[bb] += m*(dx*dvy - dy*dvx)/rr;
	dvr[bb] += m*(dx*dvx + dy*dvy)/rr;
      }
      continue;
    }
    if(nn < species[SPECIES_GAS].offset) {
      continue;
    }
    // APIC spin angular momentum of the particle's affine velocity field
    lz += m*(aff[4*(nn - species[SPECIES_GAS].offset) + 2] - aff[4*(nn - species[SPECIES_GAS].offset) + 1]);
    dx = r[nn*NUM_DIMS + 0] - cx;
    dy = r[nn*NUM_DIMS + 1] - cy;
    dvx = v[nn*NUM_DIMS + 0] - svx;
    dvy = v[nn*NUM_DIMS + 1] - svy;
    rr = sqrt(dx*dx + dy*dy);
    bb = (int)(rr);
    if(rr > 0.0 && bb < NB) {
      mass[bb] += m;
      vphi[bb] += m*(dx*dvy - dy*dvx)/rr;
      vr[bb] += m*(dx*dvx + dy*dvy)/rr;
      vr2[bb] += m*(dx*dvx + dy*dvy)*(dx*dvx + dy*dvy)/(rr*rr);
      csb[bb] += m*sqrt(GAMMA*(GAMMA - 1.0)*u[nn - species[SPECIES_GAS].offset]);
    }
    if(!sink[nn - species[SPECIES_GAS].offset]) {
      e_th += m*u[nn - species[SPECIES_GAS].offset];
    }
    else {
      m_acc += m;
    }
  }

  // v_K is the star's circular speed at R_out, cs the initial sound speed
  // in the middle of the disk
  fprintf(fp, "# t %.6f star_x %.6f star_y %.6f star_vx %.6e star_vy %.6e cm_x %.6f cm_y %.6f Lz %.10e "
	  "R_in %.3f R_out %.3f edge_out %.3f v_K %.6f cs %.6f disk_to_star %.4f G %.6f gamma %.6f beta %.4f "
	  "cooling %d Q_init %.4f Q_floor %.4f E_th %.10e E_rad %.10e M_acc %.6e",
	  t, cx, cy, svx, svy, gx, gy, lz, SD_R_IN, SD_R_OUT, SD_EDGE_OUT, SD_V_K, sd_cs_init(0.5*(SD_R_IN + SD_R_OUT)),
	  SD_DISK_TO_STAR, G_DISK, GAMMA, SD_BETA, SD_COOLING, SD_Q_INIT, SD_Q_FLOOR, e_th, e_rad, m_acc);
  if(NUM_DUST_PARTICLES > 0) {
    fprintf(fp, " dust_to_gas %.6f St %.6f M_acc_dust %.6e", (double)(DUST_TO_GAS), (double)(DUST_ST), m_acc_dust);
  }
  fprintf(fp, "\n");
  for(bb=0; bb<NB; bb++) {
    if(mass[bb] > 0.0) {
      fprintf(fp, "%.2f %.8e %.6e %.6e %.6e %.6e", bb + 0.5, mass[bb]/(M_PI*(double)(2*bb + 1)),
	      vphi[bb]/mass[bb], vr[bb]/mass[bb], sqrt(vr2[bb]/mass[bb]), csb[bb]/mass[bb]);
      if(NUM_DUST_PARTICLES > 0) {
	fprintf(fp, " %.8e %.6e %.6e", dmass[bb]/(M_PI*(double)(2*bb + 1)),
		dmass[bb] > 0.0 ? dvphi[bb]/dmass[bb] : 0.0, dmass[bb] > 0.0 ? dvr[bb]/dmass[bb] : 0.0);
      }
      fprintf(fp, "\n");
    }
  }
  fprintf(fp, "\n\n");
  fflush(fp);

  printf("disk t: %.3f star: (%.4f, %.4f) v (%.3e, %.3e) Lz %.8e E_th %.6e E_rad %.6e M_acc %.4e\n",
	 t, cx, cy, svx, svy, lz, e_th, e_rad, m_acc);
}
#endif

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
  const double *affx_gas = &grid_fields[FIELD_GAS_AFF_X*NUM_CELLS*NUM_CELLS];
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
      mx += mx_gas[nn] + affx_gas[nn];
      e += e_gas[nn];
    }
    // mesh nodes sit at integer coordinates
    fprintf(fp, "%d %.10g %.10g %.10g %.10g\n", ii,
	    rho/NUM_CELLS/ST_RHO_R, mx/rho, (GAMMA - 1.0)*e/NUM_CELLS/p_r, e/rho);
  }

  fclose(fp);
}
#endif

// dust drag fields on the mesh, NUM_DRAG_FIELDS arrays of NUM_CELLS^2
enum drag_field_id{
  DRAG_DVG_X = 0,  // gas velocity kick of the drag
  DRAG_DVG_Y,
  DRAG_DVD_X,      // dust velocity kick of the drag
  DRAG_DVD_Y,
  DRAG_DAMP,       // relaxation 1 - exp(-dt/t_s) of the dust velocity dispersion
  DRAG_VD_X,       // mass weighted dust velocity including the affine momentum
  DRAG_VD_Y,
  DRAG_HEAT,       // kinetic energy change of the drag kick, CIC deposited
  NUM_DRAG_FIELDS
};

// dust stopping time at mesh node nn with gas surface density rho_g
static double dust_stop_time(int nn, double rho_g, const double *r){
  if(DUST_TS > 0.0) {
    return DUST_TS;
  }
#if STAR_DISK
  // Epstein drag in a thin disk, St = Omega t_s ~ 1/sigma_gas
  int ns = species[SPECIES_DM].offset;
  double dx = (double)(nn % NUM_CELLS) - r[ns*NUM_DIMS + 0];
  double dy = (double)(nn / NUM_CELLS) - r[ns*NUM_DIMS + 1];

  return DUST_ST*sd_sigma_ref/(rho_g*sd_lookup(sd_omega_tab, sqrt(dx*dx + dy*dy)));
#else
  (void)nn;
  (void)rho_g;
  (void)r;
  return 1.0e30;
#endif
}

// two-fluid drag at the mesh nodes over the step dt. With dust to gas ratio
// eps = rho_d/rho_g the gas and dust velocities obey
//   dv_g/dt = a_g - eps (v_g - v_d)/t_s,   dv_d/dt = a_d - (v_d - v_g)/t_s
// The center of mass velocity only feels the forces, the relative velocity
// w = v_g - v_d relaxes on tau = t_s/(1 + eps) towards the terminal value
// set by the differential acceleration da = a_g - a_d, the gas pressure
// acceleration (gravity acts on both alike):
//   w(dt) = da tau + (w0 - da tau) exp(-dt/tau)
// The kicks of the step already apply da dt, the drag adds
//   dw = w(dt) - w0 - da dt
// shared between the species so that momentum is conserved,
// rho_g dv_g + rho_d dv_d = 0. The solution is exact for any dt/t_s, and
// keeping da in it gives the terminal drift w = da tau of tightly coupled
// grains, which a split exponential drag after the force kick would damp
// away (Lorén-Aguilar and Bate 2015, Yang and Johansen 2016). w0 is the
// difference of the mass weighted mesh velocities deposited at the half step,
// both including their APIC affine momentum
void compute_dust_drag(const double *grid_fields, const double *vel, const double *acc, double *drag,
		       double dt, const double *r){
  int nn;
  double rg, rd, ts, tau, e, wx, wy, dax, day, dwx, dwy;
  const double *rho_gas = &grid_fields[FIELD_RHO_GAS*NUM_CELLS*NUM_CELLS];
  const double *rho_dust = &grid_fields[FIELD_RHO_DUST*NUM_CELLS*NUM_CELLS];
  const double *mx_dust = &grid_fields[FIELD_DUST_MOM_X*NUM_CELLS*NUM_CELLS];
  const double *my_dust = &grid_fields[FIELD_DUST_MOM_Y*NUM_CELLS*NUM_CELLS];
  const double *affx_dust = &grid_fields[FIELD_DUST_AFF_X*NUM_CELLS*NUM_CELLS];
  const double *affy_dust = &grid_fields[FIELD_DUST_AFF_Y*NUM_CELLS*NUM_CELLS];
  const double *vgx = &vel[0];
  const double *vgy = &vel[NUM_CELLS*NUM_CELLS];
  const double *ax = &acc[0];
  const double *ay = &acc[NUM_CELLS*NUM_CELLS];
  double *f[NUM_DRAG_FIELDS];

  for(nn=0; nn<NUM_DRAG_FIELDS; nn++) {
    f[nn] = &drag[nn*NUM_CELLS*NUM_CELLS];
  }

  for(nn=0; nn<NUM_CELLS*NUM_CELLS; nn++) {
    rg = rho_gas[nn];
    rd = rho_dust[nn];
    f[DRAG_VD_X][nn] = rd > 0.0 ? (mx_dust[nn] + affx_dust[nn])/rd : 0.0;
    f[DRAG_VD_Y][nn] = rd > 0.0 ? (my_dust[nn] + affy_dust[nn])/rd : 0.0;

    if(rg <= RHO_GAS_FLOOR || rd <= 0.0 || dt <= 0.0) {
      f[DRAG_DVG_X][nn] = f[DRAG_DVG_Y][nn] = 0.0;
      f[DRAG_DVD_X][nn] = f[DRAG_DVD_Y][nn] = 0.0;
      f[DRAG_DAMP][nn] = 0.0;
      continue;
    }

    ts = dust_stop_time(nn, rg, r);
    tau = ts/(1.0 + rd/rg);
    e = exp(-dt/tau);

    wx = vgx[nn] - f[DRAG_VD_X][nn];
    wy = vgy[nn] - f[DRAG_VD_Y][nn];
    dax = ax[nn];
    day = ay[nn];
    dwx = -(1.0 - e)*(wx - dax*tau) - dax*dt;
    dwy = -(1.0 - e)*(wy - day*tau) - day*dt;

    f[DRAG_DVG_X][nn] = rd/(rg + rd)*dwx;
    f[DRAG_DVG_Y][nn] = rd/(rg + rd)*dwy;
    f[DRAG_DVD_X][nn] = -rg/(rg + rd)*dwx;
    f[DRAG_DVD_Y][nn] = -rg/(rg + rd)*dwy;
    f[DRAG_DAMP][nn] = 1.0 - exp(-dt/ts);
  }
}

// drag kick of the particles. Gas particles take the CIC interpolated gas
// kick. Dust particles take the dust kick and relax their velocity relative
// to the dust velocity of each node of their stencil,
//   dv_p = sum_c w_c (dv_d,c - f_c (v_p - V_d,c)),
// so a dust particle forgets its own velocity on its stopping time, the
// velocity dispersion of the dust that the node drag cannot see. For grains
// with t_s << dt (f = 1) this is a full PIC transfer of the dust velocity,
// whose CIC smoothing would damp resolved dust motion by ~k^2/6 per step, so
// dust carries an APIC matrix B like the gas: V_d,c includes the deposited
// affine momentum and B relaxes towards sum_c w_c V'_c (x_c - x_p)^T of the
// kicked mesh velocity V' = V_d + dv_d with the same fraction f. Summed over
// particles the relaxation adds sum_c f_c aff_c to the momentum, which vanishes
// for uniform f (the APIC terms of a particle sum to zero) and is small
// where f varies slowly. The gas kick and the dust kick conserve momentum
// exactly. Each particle deposits its kinetic energy change with its CIC
// weights into the per-thread heat subfield
void *dust_drag_thread(void *threadarg){
  int nn, lo, hi, ii, jj, kk, s;
  struct thread_data *my_data = (struct thread_data *) threadarg;
  int thread_id = my_data->thread_id;
  double *r = my_data->r, *v = my_data->v;
  const double *drag = my_data->drag_mesh;
  const double *dvgx = &drag[DRAG_DVG_X*NUM_CELLS*NUM_CELLS], *dvgy = &drag[DRAG_DVG_Y*NUM_CELLS*NUM_CELLS];
  const double *dvdx = &drag[DRAG_DVD_X*NUM_CELLS*NUM_CELLS], *dvdy = &drag[DRAG_DVD_Y*NUM_CELLS*NUM_CELLS];
  const double *damp = &drag[DRAG_DAMP*NUM_CELLS*NUM_CELLS];
  const double *vdx = &drag[DRAG_VD_X*NUM_CELLS*NUM_CELLS], *vdy = &drag[DRAG_VD_Y*NUM_CELLS*NUM_CELLS];
  double *heat = &my_data->sub_heat[thread_id*NUM_CELLS*NUM_CELLS];
  double dr[2], w[4], gw[4][2], d[4][2];
  int c[4];
  double m, vx, vy, dx, dy, dk, ke = 0.0;
  const double *rho_dust = &my_data->grid_fields[FIELD_RHO_DUST*NUM_CELLS*NUM_CELLS];
  double fp, bg[4], rho_min, sparse, *bp;

  memset(heat, 0, NUM_CELLS*NUM_CELLS*sizeof(double));

  for(s=SPECIES_GAS; s<=SPECIES_DUST; s++) {
    m = species[s].mass;
    thread_slice(species[s].offset, species[s].count, thread_id, NUM_THREADS, &lo, &hi);
    for(nn=lo; nn<hi; nn++) {
      if((s == SPECIES_GAS && my_data->sink[nn - species[SPECIES_GAS].offset]) ||
	 (s == SPECIES_DUST && my_data->dsink[nn - species[SPECIES_DUST].offset])) {
	continue;
      }
      ii = floor(r[nn*NUM_DIMS + 0]);
      jj = floor(r[nn*NUM_DIMS + 1]);
      dr[0] = r[nn*NUM_DIMS + 0] - (double)(ii);
      dr[1] = r[nn*NUM_DIMS + 1] - (double)(jj);
      cic_stencil(ii, jj, dr, NUM_CELLS, c, w, gw, d);

      vx = v[nn*NUM_DIMS + 0];
      vy = v[nn*NUM_DIMS + 1];
      dx = 0.0;
      dy = 0.0;
      if(s == SPECIES_GAS) {
	for(kk=0; kk<4; kk++) {
	  dx += w[kk]*dvgx[c[kk]];
	  dy += w[kk]*dvgy[c[kk]];
	}
      }
      else {
	fp = 0.0;
	bg[0] = bg[1] = bg[2] = bg[3] = 0.0;
	rho_min = rho_dust[c[0]];
	for(kk=0; kk<4; kk++) {
	  dx += w[kk]*(dvdx[c[kk]] - damp[c[kk]]*(vx - vdx[c[kk]]));
	  dy += w[kk]*(dvdy[c[kk]] - damp[c[kk]]*(vy - vdy[c[kk]]));
	  fp += w[kk]*damp[c[kk]];
	  bg[0] += w[kk]*(vdx[c[kk]] + dvdx[c[kk]])*d[kk][0];
	  bg[1] += w[kk]*(vdx[c[kk]] + dvdx[c[kk]])*d[kk][1];
	  bg[2] += w[kk]*(vdy[c[kk]] + dvdy[c[kk]])*d[kk][0];
	  bg[3] += w[kk]*(vdy[c[kk]] + dvdy[c[kk]])*d[kk][1];
	  rho_min = fmin(rho_min, rho_dust[c[kk]]);
	}
	// near-empty nodes get the same limiter as the gas APIC
	sparse = fmin(1.0, rho_min/(APIC_SPARSE*m));
	bp = &my_data->daff[4*(nn - species[SPECIES_DUST].offset)];
	for(kk=0; kk<4; kk++) {
	  bp[kk] += fp*(bg[kk] - bp[kk]);
	  bp[kk] *= sparse;
	}
      }
      v[nn*NUM_DIMS + 0] = vx + dx;
      v[nn*NUM_DIMS + 1] = vy + dy;

      dk = 0.5*m*(dx*(2.0*vx + dx) + dy*(2.0*vy + dy));
      ke += dk;
      for(kk=0; kk<4; kk++) {
	heat[c[kk]] += w[kk]*dk;
      }
    }
  }
  my_data->drag_ke = ke;

  pthread_exit(NULL);
}

// the kinetic energy the drag removed heats the gas: a gas particle gains
// u += -sum_c w_c H_c/rho_c from the summed heat field H. Summed over the gas
// with the deposit weights this returns sum_c H_c exactly, except at nodes
// without gas and at the internal energy floor (the residual drag_ke + drag_heat)
void *dust_heat_thread(void *threadarg){
  int nn, lo, hi, ii, jj, kk, gg;
  struct thread_data *my_data = (struct thread_data *) threadarg;
  int thread_id = my_data->thread_id;
  const double *r = my_data->r;
  double *u = my_data->u;
  const double *heat = &my_data->drag_mesh[DRAG_HEAT*NUM_CELLS*NUM_CELLS];
  const double *rho_gas = &my_data->grid_fields[FIELD_RHO_GAS*NUM_CELLS*NUM_CELLS];
  double dr[2], w[4], gw[4][2], d[4][2];
  int c[4];
  double h, u_new, eh = 0.0;

  thread_slice(species[SPECIES_GAS].offset, species[SPECIES_GAS].count, thread_id, NUM_THREADS, &lo, &hi);
  for(nn=lo; nn<hi; nn++) {
    gg = nn - species[SPECIES_GAS].offset;
    if(my_data->sink[gg]) {
      continue;
    }
    ii = floor(r[nn*NUM_DIMS + 0]);
    jj = floor(r[nn*NUM_DIMS + 1]);
    dr[0] = r[nn*NUM_DIMS + 0] - (double)(ii);
    dr[1] = r[nn*NUM_DIMS + 1] - (double)(jj);
    cic_stencil(ii, jj, dr, NUM_CELLS, c, w, gw, d);

    h = 0.0;
    for(kk=0; kk<4; kk++) {
      if(rho_gas[c[kk]] > RHO_GAS_FLOOR) {
	h -= w[kk]*heat[c[kk]]/rho_gas[c[kk]];
      }
    }
#if !ISOTHERMAL
    u_new = fmax(u[gg] + h, U_MIN);
    eh += u_new - u[gg];
    u[gg] = u_new;
#else
    (void)u_new;
    (void)u;
#endif
  }
  my_data->drag_heat = species[SPECIES_GAS].mass*eh;

  pthread_exit(NULL);
}

#if DUSTY_BOX
// displacement, velocity and adiabatic internal energy factor of the sound
// wave at lattice position x0: the density is rho0/(1 + dxi/dx0) =
// rho0 (1 + DB_AMP sin(k x)) to first order, the gas velocity of a wave
// travelling to +x is cs drho/rho and u ~ rho^(gamma - 1), so the wave
// carries no entropy perturbation
static void dusty_wave(double x0, double *x, double *vx, double *fu){
#if DB_WAVE
  *x = x0 + DB_AMP/DB_K*cos(DB_K*x0);
  *vx = DB_CS*DB_AMP*sin(DB_K*x0);
  *fu = 1.0 + (GAMMA - 1.0)*DB_AMP*sin(DB_K*x0);
#else
  *x = x0;
  *vx = 0.0;
  *fu = 1.0;
#endif
}

// uniform gas at rest and dust moving with DB_V0 along x on DB_PPC_X by
// DB_PPC_Y lattices, or both carrying the sound wave
void init_dusty_box(double *r, double *v, double *a, double *u, double *du){
  int ii, jj, px, py;
  int nn = species[SPECIES_GAS].offset;
  double x, vx, fu;

  for(jj=0; jj<NUM_CELLS; jj++) {
    for(ii=0; ii<NUM_CELLS; ii++) {
      for(py=0; py<DB_PPC_Y; py++) {
	for(px=0; px<DB_PPC_X; px++) {
	  dusty_wave(ii + (px + 0.5)/DB_PPC_X, &x, &vx, &fu);
	  r[nn*NUM_DIMS + 0] = wrap_periodic(x, (double)(NUM_CELLS));
	  r[nn*NUM_DIMS + 1] = jj + (py + 0.5)/DB_PPC_Y;
	  v[nn*NUM_DIMS + 0] = vx;
	  v[nn*NUM_DIMS + 1] = 0.0;
	  u[nn - species[SPECIES_GAS].offset] = fu*DB_CS*DB_CS/(GAMMA*(GAMMA - 1.0));
	  du[nn - species[SPECIES_GAS].offset] = 0.0;
	  nn++;
	}
      }
    }
  }
  for(jj=0; jj<NUM_CELLS; jj++) {
    for(ii=0; ii<NUM_CELLS; ii++) {
      for(py=0; py<DB_PPC_Y; py++) {
	for(px=0; px<DB_PPC_X; px++) {
	  dusty_wave(ii + (px + 0.5)/DB_PPC_X, &x, &vx, &fu);
	  r[nn*NUM_DIMS + 0] = wrap_periodic(x, (double)(NUM_CELLS));
	  r[nn*NUM_DIMS + 1] = jj + (py + 0.5)/DB_PPC_Y;
	  v[nn*NUM_DIMS + 0] = DB_V0 + vx;
	  v[nn*NUM_DIMS + 1] = 0.0;
	  nn++;
	}
      }
    }
  }
  for(nn=0; nn<NUM_PARTICLES; nn++) {
    a[nn*NUM_DIMS + 0] = 0.0;
    a[nn*NUM_DIMS + 1] = 0.0;
  }
}

// one line per step: time, mass weighted gas and dust velocity along x and
// y, total momentum along x, kinetic and thermal energy, cumulative residual
// of the drag heating
void write_dusty_box(FILE *fp, const double *v, const double *u, double t, double dt, double e_drag_res){
  int nn, s;
  double p[NUM_SPECIES][2], ke = 0.0, eth = 0.0, m;

  for(s=0; s<NUM_SPECIES; s++) {
    p[s][0] = p[s][1] = 0.0;
    m = species[s].mass;
    for(nn=species[s].offset; nn<species[s].offset + species[s].count; nn++) {
      p[s][0] += m*v[nn*NUM_DIMS + 0];
      p[s][1] += m*v[nn*NUM_DIMS + 1];
      ke += 0.5*m*(v[nn*NUM_DIMS + 0]*v[nn*NUM_DIMS + 0] + v[nn*NUM_DIMS + 1]*v[nn*NUM_DIMS + 1]);
    }
  }
  for(nn=0; nn<NUM_GAS_PARTICLES; nn++) {
    eth += GAS_MASS*u[nn];
  }
  if(t == 0.0) {
    fprintf(fp, "# eps %.10g t_s %.10g v0 %.10g dt %.10g gamma %.10g\n", DB_EPS, DUST_TS, DB_V0, dt, GAMMA);
    fprintf(fp, "# t v_gas_x v_gas_y v_dust_x v_dust_y p_x E_kin E_th drag_res\n");
  }
  fprintf(fp, "%.10e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.6e\n", t,
	  p[SPECIES_GAS][0]/(GAS_MASS*NUM_GAS_PARTICLES), p[SPECIES_GAS][1]/(GAS_MASS*NUM_GAS_PARTICLES),
	  p[SPECIES_DUST][0]/(DUST_MASS*NUM_DUST_PARTICLES), p[SPECIES_DUST][1]/(DUST_MASS*NUM_DUST_PARTICLES),
	  p[SPECIES_GAS][0] + p[SPECIES_DUST][0], ke, eth, e_drag_res);
  fflush(fp);
}

#if DB_WAVE
// one line per step: time and the Fourier coefficients at the wave number
// DB_K of the gas and dust density and momentum, normalized by the species
// mass, sum_p m exp(-i k x_p)/M and sum_p m v_p exp(-i k x_p)/M as real and
// imaginary parts
void write_dusty_wave(FILE *fp, const double *r, const double *v, const double *u, double t){
  int nn, s, kk = 0;
  double out[8], ph, m, mt, eth = 0.0;

  for(s=SPECIES_GAS; s<=SPECIES_DUST; s++) {
    double dr = 0.0, di = 0.0, pr = 0.0, pi = 0.0;
    m = species[s].mass;
    mt = m*(double)(species[s].count);
    for(nn=species[s].offset; nn<species[s].offset + species[s].count; nn++) {
      ph = DB_K*r[nn*NUM_DIMS + 0];
      dr += cos(ph);
      di -= sin(ph);
      pr += v[nn*NUM_DIMS + 0]*cos(ph);
      pi -= v[nn*NUM_DIMS + 0]*sin(ph);
    }
    out[kk++] = m*dr/mt;
    out[kk++] = m*di/mt;
    out[kk++] = m*pr/mt;
    out[kk++] = m*pi/mt;
  }
  for(nn=0; nn<NUM_GAS_PARTICLES; nn++) {
    eth += GAS_MASS*u[nn];
  }
  if(t == 0.0) {
    fprintf(fp, "# eps %.10g t_s %.10g cs %.10g k %.10g amp %.10g gamma %.10g dt %.10g\n",
	    DB_EPS, DUST_TS, DB_CS, DB_K, DB_AMP, GAMMA, DT);
    fprintf(fp, "# t gas_rho_re gas_rho_im gas_mom_re gas_mom_im dust_rho_re dust_rho_im dust_mom_re dust_mom_im E_th\n");
  }
  fprintf(fp, "%.10e", t);
  for(kk=0; kk<8; kk++) {
    fprintf(fp, " %.12e", out[kk]);
  }
  fprintf(fp, " %.15e\n", eth);
  fflush(fp);
}
#endif
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

  
  double *rho, *sub_rho, *rho_species;
  double *u, *du, *aff, *daff;
  unsigned char *sink, *dsink;
  double *press, *du_mesh, *gas_vel, *gas_acc;
  double *drag_mesh, *sub_heat;
  double e_drag_ke = 0.0, e_drag_res = 0.0;
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
  aff=(double *)calloc(4*(NUM_GAS_PARTICLES > 0 ? NUM_GAS_PARTICLES : 1), sizeof(double));
  sink=(unsigned char *)calloc(NUM_GAS_PARTICLES > 0 ? NUM_GAS_PARTICLES : 1, sizeof(unsigned char));
  dsink=(unsigned char *)calloc(NUM_DUST_PARTICLES > 0 ? NUM_DUST_PARTICLES : 1, sizeof(unsigned char));
  daff=(double *)calloc(4*(NUM_DUST_PARTICLES > 0 ? NUM_DUST_PARTICLES : 1), sizeof(double));
  if(!u || !du || !aff || !sink || !dsink || !daff){
    printf("Out of memory: u or du not allocated.\n");
    exit(1);
  }
  // gas pressure, du/dt and mesh velocity
  press=(double *)malloc(NUM_CELLS*NUM_CELLS*sizeof(double));
  du_mesh=(double *)malloc(NUM_CELLS*NUM_CELLS*sizeof(double));
  gas_vel=(double *)malloc(4*NUM_CELLS*NUM_CELLS*sizeof(double));
  gas_acc=(double *)malloc(2*NUM_CELLS*NUM_CELLS*sizeof(double));
  if(!press || !du_mesh || !gas_vel || !gas_acc){
    printf("Out of memory: gas pressure fields not allocated.\n");
    exit(1);
  }
  // dust drag fields and per-thread heat subfields
#if NUM_DUST_PARTICLES > 0
  drag_mesh=(double *)malloc(NUM_DRAG_FIELDS*NUM_CELLS*NUM_CELLS*sizeof(double));
  sub_heat=(double *)malloc(NUM_THREADS*NUM_CELLS*NUM_CELLS*sizeof(double));
  if(!drag_mesh || !sub_heat){
    printf("Out of memory: dust drag fields not allocated.\n");
    exit(1);
  }
#else
  drag_mesh=NULL;
  sub_heat=NULL;
  (void)drag_mesh;
  (void)sub_heat;
  (void)e_drag_ke;
  (void)e_drag_res;
#endif
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
#elif DUSTY_BOX
  init_dusty_box(r, v, a, u, du);
#else
  init_star_disk(r, v, a, u, du);
#endif
  cs_max = 0.0;
  v_max = 0.0;
  
  // compute gravitational green's function in k-space for the 5 point
  // discrete laplacian, FFT bin ii holds wavenumber ii for ii <= N/2 and
  // ii - N above, the k = 0 mean density mode is dropped (periodic box).
  // The solve multiplies by 4 pi G_EFF, so -1/K^2 with K^2 = 4 dd_k gives
  // lap(phi) = 4 pi G_EFF rho. Thin disk gravity is solved separately below
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

#if THIN_DISK_GRAVITY
  // isolated thin disk gravity (Hockney and Eastwood): the density is zero
  // padded to a mesh of twice the size and convolved with the Plummer
  // softened point mass potential -G_DISK/sqrt(r^2 + h^2), the transform of
  // phi_k = -2 pi G_DISK exp(-k h) sigma_k/k, tabulated with the distances
  // wrapped on the padded mesh. No periodic images act on the disk
  fftw_complex *pad_in, *pad_hat, *pad_out, *kern_hat;
  fftw_plan pad_fwd, pad_bwd;
  const int NP = 2*NUM_CELLS;
  double pdx, pdy;

  pad_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NP*NP);
  pad_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NP*NP);
  pad_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NP*NP);
  kern_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NP*NP);
  if(!pad_in || !pad_hat || !pad_out || !kern_hat){
    printf("Out of memory: padded gravity mesh not allocated.\n");
    exit(1);
  }
  pad_fwd = fftw_plan_dft_2d(NP, NP, pad_in, pad_hat, FFTW_FORWARD, FFTW_ESTIMATE);
  pad_bwd = fftw_plan_dft_2d(NP, NP, pad_hat, pad_out, FFTW_BACKWARD, FFTW_ESTIMATE);
  if(!pad_fwd || !pad_bwd){
    printf("FFTW padded gravity plans failed.\n");
    exit(1);
  }
  for(jj=0; jj<NP; jj++) {
    for(ii=0; ii<NP; ii++) {
      pdx = (double)(ii <= NP/2 ? ii : NP - ii);
      pdy = (double)(jj <= NP/2 ? jj : NP - jj);
      pad_in[jj*NP + ii][0] = -G_DISK/sqrt(pdx*pdx + pdy*pdy + SD_SOFT*SD_SOFT);
      pad_in[jj*NP + ii][1] = 0.0;
    }
  }
  fftw_execute(pad_fwd);
  memcpy(kern_hat, pad_hat, sizeof(fftw_complex)*NP*NP);
#endif

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

#if STAR_DISK
  // the first pass through the loop runs with dt = 0 to compute the forces
  // of the initial positions, then the rotation is set from them
  int setup_pass = 1;
  double next_disk_out = SD_OUT_DT;
  double e_rad = 0.0;
  FILE *disk_fp = fopen("stardisk.dat", "w");
  FILE *map_fp = fopen("stardisk_rho.bin", "wb");
  FILE *dust_map_fp = NUM_DUST_PARTICLES > 0 ? fopen("stardisk_dust.bin", "wb") : NULL;
  if(!disk_fp || !map_fp || (NUM_DUST_PARTICLES > 0 && !dust_map_fp)) {
    printf("Could not open stardisk.dat, stardisk_rho.bin or stardisk_dust.bin for writing.\n");
    exit(1);
  }
  printf("star disk: disk/star mass %.3f, R %.1f - %.1f cells, Q_init %.2f, cooling %d beta %.2f Q_floor %.2f\n"
	 "           cs %.4f and orbit %.4f at R_in, cs %.4f and orbit %.4f at R_out\n",
	 SD_DISK_TO_STAR, SD_R_IN, SD_R_OUT, SD_Q_INIT, SD_COOLING, SD_BETA, SD_Q_FLOOR,
	 sd_cs_init(SD_R_IN + SD_EDGE_IN), 2.0*M_PI/sd_lookup(sd_omega_tab, SD_R_IN),
	 sd_cs_init(SD_R_OUT), 2.0*M_PI/sd_lookup(sd_omega_tab, SD_R_OUT));
  if(NUM_DUST_PARTICLES > 0) {
    printf("           dust: %d particles, dust/gas mass %.4f, St %.3f at the initial gas surface density at r %.1f\n",
	   NUM_DUST_PARTICLES, (double)(DUST_TO_GAS), (double)(DUST_ST), 0.5*(SD_R_IN + SD_R_OUT));
  }
  dt = 0.0;
#endif
  
#if DUSTY_BOX
  FILE *box_fp = fopen(DB_WAVE ? "dustywave.dat" : "dustybox.dat", "w");
  if(!box_fp) {
    printf("Could not open the dusty box output for writing.\n");
    exit(1);
  }
#if DB_WAVE
  write_dusty_wave(box_fp, r, v, u, 0.0);
#else
  write_dusty_box(box_fp, v, u, 0.0, dt, 0.0);
#endif
#endif

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
      thread_data_array[nn].aff=aff;
      thread_data_array[nn].sink=sink;
      thread_data_array[nn].dsink=dsink;
      thread_data_array[nn].daff=daff;
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
    
#if THIN_DISK_GRAVITY
    // isolated solve: zero padded density convolved with the softened
    // point mass potential, rho holds the mass per cell
    memset(pad_in, 0, sizeof(fftw_complex)*NP*NP);
    for(jj=0; jj<NUM_CELLS; jj++) {
      for(ii=0; ii<NUM_CELLS; ii++) {
	pad_in[jj*NP + ii][0] = rho[jj*NUM_CELLS + ii];
      }
    }
    fftw_execute(pad_fwd);
    for(ii=0; ii<NP*NP; ii++) {
      double re = pad_hat[ii][0]*kern_hat[ii][0] - pad_hat[ii][1]*kern_hat[ii][1];
      double im = pad_hat[ii][0]*kern_hat[ii][1] + pad_hat[ii][1]*kern_hat[ii][0];
      pad_hat[ii][0] = re;
      pad_hat[ii][1] = im;
    }
    fftw_execute(pad_bwd);
    for(jj=0; jj<NUM_CELLS; jj++) {
      for(ii=0; ii<NUM_CELLS; ii++) {
	phi[jj*NUM_CELLS + ii] = pad_out[jj*NP + ii][0]/((double)(NP)*(double)(NP));
      }
    }
#else
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
#endif

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
      thread_data_array[nn].aff=aff;
      thread_data_array[nn].sink=sink;
    }
    run_threads(hydro_interpolation_thread);
    for(nn=0;nn<NUM_THREADS;nn++){
      e_flip_loss += thread_data_array[nn].flip_loss;
    }
    e_work_res += p_res*dt;

#if NUM_DUST_PARTICLES > 0
    // gas-dust drag: exact two-fluid solution at the mesh nodes, kicks of the
    // particles and the dissipated kinetic energy returned to the gas
    compute_dust_drag(rho_species, gas_vel, gas_acc, drag_mesh, dt, r);
    for(nn=0;nn<NUM_THREADS;nn++){
      thread_data_array[nn].thread_id=nn;
      thread_data_array[nn].r=r;
      thread_data_array[nn].v=v;
      thread_data_array[nn].u=u;
      thread_data_array[nn].sink=sink;
      thread_data_array[nn].dsink=dsink;
      thread_data_array[nn].daff=daff;
      thread_data_array[nn].drag_mesh=drag_mesh;
      thread_data_array[nn].sub_heat=sub_heat;
      thread_data_array[nn].grid_fields=rho_species;
    }
    run_threads(dust_drag_thread);
    for(ii=0;ii<NUM_CELLS*NUM_CELLS;ii++){
      acc = 0.0;
      for(nn=0;nn<NUM_THREADS;nn++){
	acc += sub_heat[nn*NUM_CELLS*NUM_CELLS + ii];
      }
      drag_mesh[DRAG_HEAT*NUM_CELLS*NUM_CELLS + ii] = acc;
    }
    run_threads(dust_heat_thread);
    for(nn=0;nn<NUM_THREADS;nn++){
      e_drag_ke += thread_data_array[nn].drag_ke;
      e_drag_res += thread_data_array[nn].drag_ke + thread_data_array[nn].drag_heat;
    }
#endif

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

#if STAR_DISK
    if(setup_pass) {
      v_max = set_equilibrium_rotation(r, v, a);
      setup_pass = 0;
      dt = fmin(DT, CFL_TARGET/(cs_max + v_max));
      write_star_disk_profile(disk_fp, r, v, aff, u, sink, dsink, e_rad, 0.0);
      write_star_disk_map(map_fp, &rho_species[FIELD_RHO_GAS*NUM_CELLS*NUM_CELLS]);
      if(dust_map_fp) {
	write_star_disk_map(dust_map_fp, &rho_species[FIELD_RHO_DUST*NUM_CELLS*NUM_CELLS]);
      }
      continue;
    }
#endif

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
	  // dust shares the blue channel with the collisionless species
	  d_dm=display_scale[SPECIES_DM]*rho_species[SPECIES_DM*NUM_CELLS*NUM_CELLS + jj*NUM_CELLS+ii] +
	    display_scale[SPECIES_DUST]*rho_species[SPECIES_DUST*NUM_CELLS*NUM_CELLS + jj*NUM_CELLS+ii];
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
#if NUM_DUST_PARTICLES > 0
      printf("dust drag: kinetic energy change %.6e, not returned as heat %.3e\n", e_drag_ke, e_drag_res);
#endif

#if TIFF_ENABLE
      sprintf(filename, "/home/janne808/testrun/%08d.tif", tiff_frame++);
      writeframe(filename, screen);
#endif

      // timer start
      clock_gettime(CLOCK_MONOTONIC, &time1);

      // reset steps variable
      steps = -1;
    }
#if STAR_DISK
#if SD_COOLING
    // radiative cooling over the completed step
    e_rad += cool_star_disk(r, u, sink, dt);
#endif
    // accretion of gas onto the star
    accrete_star_disk(r, v, aff, sink, dsink, daff);
#endif

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

#if STAR_DISK
    if(sim_time >= next_disk_out - 1.0e-12) {
      write_star_disk_profile(disk_fp, r, v, aff, u, sink, dsink, e_rad, sim_time);
      write_star_disk_map(map_fp, &rho_species[FIELD_RHO_GAS*NUM_CELLS*NUM_CELLS]);
      if(dust_map_fp) {
	write_star_disk_map(dust_map_fp, &rho_species[FIELD_RHO_DUST*NUM_CELLS*NUM_CELLS]);
      }
      next_disk_out += SD_OUT_DT;
    }
    if(sim_time + dt > next_disk_out) {
      dt = next_disk_out - sim_time;
    }
#endif

#if DUSTY_BOX
#if DB_WAVE
    write_dusty_wave(box_fp, r, v, u, sim_time);
#else
    write_dusty_box(box_fp, v, u, sim_time, dt, e_drag_res);
#endif
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
#if DUSTY_BOX
  fclose(box_fp);
#endif
#if STAR_DISK
  fclose(disk_fp);
  fclose(map_fp);
  if(dust_map_fp) {
    fclose(dust_map_fp);
  }
#endif

  // clean up FFTW
#if THIN_DISK_GRAVITY
  fftw_destroy_plan(pad_fwd);
  fftw_destroy_plan(pad_bwd);
  fftw_free(pad_in);
  fftw_free(pad_hat);
  fftw_free(pad_out);
  fftw_free(kern_hat);
#endif
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

