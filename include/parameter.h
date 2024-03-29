#ifndef SPH_PARAMETER_H
#define SPH_PARAMETER_H

//only 3 available
#define DIM 3

#define KERNEL_DIM 2

#define NN 1 // decide if nearest neighbors of each particle are also written to file

/* 0 calculate density via Kernel sum
   1 calculate density via continuity equation, then it is necessary to provide an inital density distribution */
#define CALC_DENSITY 1     

// decide whether timestep is limited by CFL number
#define COURANT_CONDITION 1

/*  0 for Leapfrog (only integrates v and r, rho calc by Kernelsum)
   1 for predictor-corrector integration 
#define INTEGRATOR 1 */

/* 0 Heun like normal scheme
   1 Heun like in Elastics Paper (no predictor-correction for r)*/
#define HEUN_LIKE_PAPER 0

/* for simulating Solids set to 1, if 0 assuming gas*/
#define SOLIDS 1

/* using artificial stress, only po*/
#define ARTIFICIAL_STRESS 1

/* use artificial viscosity or not */
#define ARTIFICIAL_VISCOSITY 1

#define DEBUG_LEVEL 1

#endif // SPH_PARAMETER_H