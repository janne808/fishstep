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

// number of particles
#define NUM_PARTICLES 2*512000

// minimum and maximum geometries
#define XMIN 0
#define XMAX NUM_CELLS

// number of dimensions
#define NUM_DIMS 2

// particle mass
#define MASS (double)(1.0/(double)(NUM_PARTICLES));

#define DT 0.05;

#define G 4.0*M_PI*M_PI;

struct thread_data thread_data_array[NUM_THREADS];

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

// vector norm
double vector_norm(double *r, int dim) {
  int ii;
  double s=0;
  
  for(ii=0; ii<dim; ii++) {
    s += r[ii]*r[ii];
  }

  return sqrt(s);
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
    if(r[ii*NUM_DIMS + 0] > (double)(NUM_CELLS)) {
      r[ii*NUM_DIMS + 0] -= (double)(NUM_CELLS-1);
    }
    else if(r[ii*NUM_DIMS + 0] < (double)(1.0)) {
      r[ii*NUM_DIMS + 0] += (double)(NUM_CELLS-1);
    }
    if(r[ii*NUM_DIMS + 1] > (double)(NUM_CELLS)) {
      r[ii*NUM_DIMS + 1] -= (double)(NUM_CELLS-1);
    }
    else if(r[ii*NUM_DIMS + 1] < (double)(1.0)) {
      r[ii*NUM_DIMS + 1] += (double)(NUM_CELLS-1);
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

    // cloud in cell field interpolation
    dr[0] = r[nn*NUM_DIMS + 0] - (double)(ii) - 0.5;
    dr[1] = r[nn*NUM_DIMS + 1] - (double)(jj) - 0.5;
    
    rho[jj*NUM_CELLS + ii] += (1.0 - dr[0])*(1.0 - dr[1])*MASS;
    rho[jj*NUM_CELLS + ((ii+1) % NUM_CELLS)] += (dr[0])*(1.0 - dr[1])*MASS;
    rho[((jj+1) % NUM_CELLS)*NUM_CELLS + ii] += (1.0 - dr[0])*(dr[1])*MASS;
    rho[((jj+1) % NUM_CELLS)*NUM_CELLS + ((ii+1) % NUM_CELLS)] += (dr[0])*(dr[1])*MASS; 
  }
}

void compute_acceleration(double *r, double *a, double *phi) {
  int nn;
  int ii, jj;

  double dr[2];
  double g_x, g_y;

  for(nn=0; nn<NUM_PARTICLES; nn++) {
    ii = floor(r[nn*NUM_DIMS + 0]);
    jj = floor(r[nn*NUM_DIMS + 1]);
    dr[0] = r[nn*NUM_DIMS + 0] - (double)(ii) - 0.5;
    dr[1] = r[nn*NUM_DIMS + 1] - (double)(jj) - 0.5;

    // calculate force grid interpolation for 2nd order differences
    g_x = -1.0*(phi[jj*NUM_CELLS + ((ii-1)%NUM_CELLS)] - phi[jj*NUM_CELLS + ((ii+1)%NUM_CELLS)])/2.0 * (1.0 - dr[0])*(1.0 - dr[1]) -
                  (phi[jj*NUM_CELLS + ii] - phi[jj*NUM_CELLS + ((ii+2)%NUM_CELLS)])/2.0 * (dr[0])*(1.0-dr[1]) -
                     (phi[((jj+1)%NUM_CELLS)*NUM_CELLS + ((ii-1)%NUM_CELLS)] - phi[((jj+1)%NUM_CELLS)*NUM_CELLS + ((ii+1)%NUM_CELLS)])/2.0 * (1.0 - dr[0])*(dr[1]) -
                        (phi[((jj+1)%NUM_CELLS)*NUM_CELLS + ii] - phi[((jj+1)%NUM_CELLS)*NUM_CELLS + ((ii+2)%NUM_CELLS)])/2.0 * (dr[0])*(dr[1]);

    g_y = -1.0*(phi[((jj-1)%NUM_CELLS)*NUM_CELLS + ii] - phi[((jj+1)%NUM_CELLS)*NUM_CELLS + ii])/2.0 * (1.0 - dr[0])*(1.0 - dr[1]) -
                  (phi[((jj-1)%NUM_CELLS)*NUM_CELLS + ((ii+1)%NUM_CELLS)] - phi[((jj+1)%NUM_CELLS)*NUM_CELLS + ((ii+1)%NUM_CELLS)])/2.0 * (dr[0])*(1.0-dr[1]) -
                     (phi[jj*NUM_CELLS + ii] - phi[((jj+2)%NUM_CELLS)*NUM_CELLS + ii])/2.0 * (1.0 - dr[0])*(dr[1]) -
                        (phi[jj*NUM_CELLS + ((ii+1)%NUM_CELLS)] - phi[((jj+2)%NUM_CELLS)*NUM_CELLS + ((ii+1)%NUM_CELLS)])/2.0 * (dr[0])*(dr[1]);
    
    a[nn*NUM_DIMS + 0] = -g_x;
    a[nn*NUM_DIMS + 1] = -g_y;
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

  // initial condition variables
  double r2[NUM_DIMS];
  double dd;
  double vv;
  
  double *rho;
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
  for(nn=0; nn<NUM_PARTICLES; nn++) {
    r[nn*NUM_DIMS + 0] = NUM_CELLS*(double)(rand())/RAND_MAX;
    r[nn*NUM_DIMS + 1] = NUM_CELLS*(double)(rand())/RAND_MAX;
    r2[0] = r[nn*NUM_DIMS + 0] - NUM_CELLS/2;
    r2[1] = r[nn*NUM_DIMS + 1] - NUM_CELLS/2;
    while(vector_norm(r2, 2) > 50.0){
      r[nn*NUM_DIMS + 0] = NUM_CELLS*(double)(rand())/RAND_MAX;
      r[nn*NUM_DIMS + 1] = NUM_CELLS*(double)(rand())/RAND_MAX;
      r2[0] = r[nn*NUM_DIMS + 0] - NUM_CELLS/2;
      r2[1] = r[nn*NUM_DIMS + 1] - NUM_CELLS/2;
    }

    dd = vector_norm(r2, 2);
    vv = sqrt(2.0*4.0*M_PI*M_PI*0.125*dd);
    
    v[nn*NUM_DIMS + 0] = -r2[1]*vv/dd;
    v[nn*NUM_DIMS + 1] = r2[0]*vv/dd;

    a[nn*NUM_DIMS + 0] = 0.0;
    a[nn*NUM_DIMS + 1] = 0.0;
  }

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

    // copy rho to complex array
    for(ii=0; ii<NUM_CELLS*NUM_CELLS; ii++) {
      rho_complex[ii][0] = (double)(4.0*M_PI*M_PI)*rho[ii];
      rho_complex[ii][1] = 0.0;
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

    compute_acceleration(r, a, phi);
    
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
	  d=NUM_PARTICLES*4.0*rho[jj*NUM_CELLS+ii];
	  //d=-0.2*phi[jj*NUM_CELLS+ii];

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
	      //colormap(d, col);
	      //pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+0]=(Uint8)(col->r*255.0);
	      //pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+1]=(Uint8)(col->g*255.0);
	      //pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+2]=(Uint8)(col->b*255.0);
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+0]=(Uint8)(d);
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+1]=(Uint8)(d);
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+2]=(Uint8)(d);
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

  fftw_destroy_plan(rho_plan);
  
#if SDL
  // clean up SDL
  SDL_Quit();
#endif

  return 0;
}

