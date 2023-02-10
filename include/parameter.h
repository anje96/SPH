#ifndef SPH_PARAMETER_H
#define SPH_PARAMETER_H

//only 3 available
#define DIM 3

#define KERNEL_DIM 2

/* 0 calculate density via Kernel sum
   1 calculate density via continuity equation, then it is necessary to provide an inital density distribution */
#define CALC_DENSITY 1            

/* 0 for Leapfrog (only integrates v and r, rho calc by Kernelsum)
   1 for predictor-corrector integration */
#define INTEGRATOR 1

/* 0 Heun like normal scheme
   1 Heun like in Elastics Paper (no predictor-correction for r)*/
#define HEUN_LIKE_PAPER 0

/* for simulating Solids set to 1, if 0 assuming gas*/
#define SOLIDS 1

#endif // SPH_PARAMETER_H