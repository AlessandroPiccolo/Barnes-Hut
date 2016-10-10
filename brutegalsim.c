/*
 * File: galsim.c
 * --------------
 * "Brute force" simulation of a galaxy
 *
 * p[0, 1, 2, 3, 4] = [mass, pos x, pos y, velocity x, velocity y]
 *
 * HPC Assignment 1
 * Alessandro Piccolo
 * Anton Sjöberg
 * Carl Christian Kirchman
 */

#include <stdio.h>
#include <stdlib.h>
#include "graphics.h"
#include "math.h"
#include "file_operations.h"
#include <sys/time.h>

const double circleRadius = 0.0015, circleColor = 0;
const int windowWidth = 1200;

static double get_wall_seconds();

int main(int argc, char *argv[]) {
  if (argc != 6) { /* Five input arguments check */
    printf("Should be,\n./galsim N input_data/filename");
    printf(" nsteps delta_t theta_max\n");
    return 1;
  }
 
  int i, j, k;                              /* Looping variables              */
  int    const N = atoi(argv[1]);           /* Number of particles            */
  int    const nsteps = atoi(argv[3]);      /* Number of steps                */
  double const dt = atof(argv[4]);          /* Time step delta t              */
  double const theta_max = atof(argv[5]);   /* Theta max threshold            */
  double *p;                                /* Buffer with particle info      */
  double  Fx[N], Fy[N];                     /* Calculated force on a particle */    
  double  r_bold[2] ;                       /* Vector between two particles   */ 
  double  r;                                /* Distance between two particles */
  double a[2];                              /* Acceleration of particle       */        
  double const G = 100/(double)N;           /* Gravitational constant         */
  double const L = 1, W = 1, epsilon = 1e-3;
  double mj,mi;

  p = (double *)malloc(N*5*sizeof(double));
  read_doubles_from_file(N*5, p, argv[2]);  /* Read particle info into buffer */
  /* Print */
  /*
  InitializeGraphics(argv[0], windowWidth, windowWidth);
  SetCAxes(0, 1);
  */
  double startTime = get_wall_seconds();
  for (k = 0; k < nsteps; k++) { /* Time steps */
    /* Calculate the force */
    for(i = 0; i < N; i++) {     /* Force on ith star from N-1 stars */
      mi = p[i*5];
      Fx[i]=0.0;
      Fy[i]=0.0;
      for(j = 0; j < N; j++) {
        if (i != j) {            /* Particle should not compare with itself */
          mj = p[j*5];
          r_bold[0] = p[i*5+1] - p[j*5+1];
          r_bold[1] = p[i*5+2] - p[j*5+2];
          r         = sqrt( r_bold[0]*r_bold[0] + r_bold[1]*r_bold[1] );
          /* Always calculate F with Plummer spheres */
          Fx[i] -=  G*mi*mj*r_bold[0]/((r+epsilon)*(r+epsilon)*(r+epsilon));
          Fy[i] -=  G*mi*mj*r_bold[1]/((r+epsilon)*(r+epsilon)*(r+epsilon));
        }
      } 
      a[0] = Fx[i]/mi;              /* Acceleration x       */
      a[1] = Fy[i]/mi;              /* Acceleration y       */
      p[i*5 + 3] += dt*a[0];        /* Update velocity x    */
      p[i*5 + 4] += dt*a[1];        /* Update velocity y    */
    } // p[0, 1, 2, 3, 4] = [mass, pos x, pos y, velocity x, velocity y]
    
    /* Update coordinates and velocity of particles */
    for (i = 0; i < N; i++) {
      p[i*5 + 1] += dt*p[i*5 + 3];  /* Update coordinates x */
      p[i*5 + 2] += dt*p[i*5 + 4];  /* Update coordinates y */
    }
    /* Print */
    /*
    if (k%10 == 0) {
      ClearScreen();  
      for(i=0;i<N;i++) { // Draw circles
      DrawCircle(p[i*5+1],p[i*5+2],L,W,circleRadius,circleColor);
      }
      Refresh();
      usleep(100); // Sleep a short while to avoid screen flickering
    }
    */
  }
  double endTime = get_wall_seconds() - startTime;
  printf("Total time = %f \n", endTime);
  write_doubles_to_file(N*5, p, "Bruteoutput.gal");

  /* Print */
  /*
  FlushDisplay();
  CloseDisplay();
  */
  free(p);
  return 0;
}

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}
