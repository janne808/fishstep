/*
 *  (C) 2019 Janne Heikkarainen <janne808@radiofreerobotron.net>
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
#define FRAMESKIP 16

// graphics scale up factor
#define SCALE 2

// number of cells per dimension
#define NUM_CELLS 256

// number of particles
#define NUM_PARTICLES 256

// minimum and maximum geometries
#define XMIN 0
#define XMAX NUM_CELLS

// number of dimensions
#define NUM_DIMS 2

// particle mass
#define MASS (double)(1.0/(double)(NUM_PARTICLES));

#define DT 0.0125;

struct thread_data thread_data_array[NUM_THREADS];

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

// particle kick-drift thread 
void *integration_kick_drift_thread(void *threadarg){
  // loop variables
  int ii;

  // thread slice boundaries
  int lo;
  int hi;

  // thread variables
  double *r, *v, *a;
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
    if(r[ii*NUM_DIMS + 0] > (double)(XMAX-1)) {
      r[ii*NUM_DIMS + 0] -= (double)(XMAX-1);
    }
    else if(r[ii*NUM_DIMS + 0] < (double)(XMIN)) {
      r[ii*NUM_DIMS + 0] += (double)(XMAX-1);
    }
    if(r[ii*NUM_DIMS + 1] > (double)(XMAX-1)) {
      r[ii*NUM_DIMS + 1] -= (double)(XMAX-1);
    }
    else if(r[ii*NUM_DIMS + 1] < (double)(XMIN)) {
      r[ii*NUM_DIMS + 1] += (double)(XMAX-1);
    }
  }

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
  double *r, *v, *a;
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

  thread_id=my_data->thread_id;

  // compute thread slice
  lo=thread_id*NUM_PARTICLES/NUM_THREADS;
  hi=(thread_id+1)*NUM_PARTICLES/NUM_THREADS;

  for(ii=lo; ii<hi; ii++) {
    v[ii*NUM_DIMS + 0] += 0.5*dt*a[ii*NUM_DIMS + 0];
    v[ii*NUM_DIMS + 1] += 0.5*dt*a[ii*NUM_DIMS + 1];
  }

  pthread_exit(NULL);
}

void compute_rho_field(double *r, double *rho) {
  // loop variables
  int nn;
  int ii, jj;

  // interpolation variables
  double dr[2];
  
  // zero rho field
  for(nn=0; nn<NUM_CELLS*NUM_CELLS; nn++) {
    rho[nn] = 0.0;
  }

  // interpolate particles into cells
  for(nn=0; nn<NUM_PARTICLES; nn++) {
    ii = floor(r[nn*NUM_DIMS + 0]);
    jj = floor(r[nn*NUM_DIMS + 1]);

    // CIC field interpolation
    dr[0] = r[nn*NUM_DIMS + 0] - (double)(ii);
    dr[1] = r[nn*NUM_DIMS + 1] - (double)(jj);
    rho[jj*NUM_CELLS + ii] += (1.0 - dr[0])*(1.0 - dr[1])*MASS;
    rho[jj*NUM_CELLS + ((ii+1) % NUM_CELLS)] += (dr[0])*(1.0 - dr[1])*MASS;
    rho[((jj+1) % NUM_CELLS)*NUM_CELLS + ii] += (1.0 - dr[0])*(dr[1])*MASS;
    rho[((jj+1) % NUM_CELLS)*NUM_CELLS + ((ii+1) % NUM_CELLS)] += (dr[0])*(dr[1])*MASS; 
  }
}

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
  int steps;
  int ii;
  int jj;
  
  // simulation variables
  double *r;
  double *v;
  double *a;

  double *rho;
  
  // time step variable
  double dt=DT;

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
  srand(1234);
  
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

  // generate initial condition
  for(nn=0; nn<NUM_PARTICLES; nn++) {
    r[nn*NUM_DIMS + 0] = NUM_CELLS*(double)(rand())/RAND_MAX;
    r[nn*NUM_DIMS + 1] = NUM_CELLS*(double)(rand())/RAND_MAX;

    v[nn*NUM_DIMS + 0] = 1.0*(2.0*(double)(rand())/RAND_MAX-1.0);
    v[nn*NUM_DIMS + 1] = 1.0*(2.0*(double)(rand())/RAND_MAX-1.0);

    a[nn*NUM_DIMS + 0] = 0.0;
    a[nn*NUM_DIMS + 1] = 0.0;
  }

  // start timer
  clock_gettime(CLOCK_MONOTONIC, &time1);

  // integrate time steps
  done=0;
  tt=0;
  
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

    // interpolate rho field from particles
    compute_rho_field(r, rho);
    
    // create kick integration threads
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(nn=0;nn<NUM_THREADS;nn++){
      thread_data_array[nn].thread_id=nn;
      thread_data_array[nn].r=r;
      thread_data_array[nn].v=v;
      thread_data_array[nn].a=a;
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
    
      for(jj=0; jj<NUM_CELLS; jj++){
	for(ii=0; ii<NUM_CELLS; ii++){
	  // get field value
	  d=NUM_PARTICLES*256.0*rho[jj*NUM_CELLS+ii];

	  // clip value
	  if(d>255.0) {
	    d=255.0;
	  }
	  else if(d<0.0) {
	    d=0.0;
	  }
	  
	  // update framebuffer
	  for(jj2=0; jj2<SCALE; jj2++){
	    for(ii2=0; ii2<SCALE; ii2++){
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+0]=(Uint8)d;
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+1]=(Uint8)d;
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+2]=(Uint8)d;
	    }
	  }
	}
      }
    
      SDL_UnlockSurface(screen);
      SDL_Flip(screen);
#endif
    
      timediff(time1, time2, &int_time);

      // print state
      printf("tt: %d steps/s: %d\n", tt, (int)(steps/((double)(int_time.tv_sec)+(double)(int_time.tv_nsec)*1.0E-9)));

#if TIFF_ENABLE
      sprintf(filename, "/home/janne808/testrun/%08d.tif", tiff_frame++);
      writeframe(filename, screen);
#endif

      // timer start
      clock_gettime(CLOCK_MONOTONIC, &time1);

      // reset steps variable
      steps = -1;
    }
    // update time step variable
    tt++;
    steps++;
  }

#if SDL
  // clean up SDL
  SDL_Quit();
#endif

  return 0;
}

