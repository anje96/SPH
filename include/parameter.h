#ifndef SPH_PARAMETER_H
#define SPH_PARAMETER_H

#define DIM 3

#define KERNEL_DIM 3

/* 0 calculate density via Kernel sum
   1 calculate density via continuity equation, then it is necessary to provide an inital density distribution */
#define CALC_DENSITY 0            

/* 0 for Leapfrog (only integrates v and r, rho calc by Kernelsum)
   1 for predictor-corrector integration (for now only integration for v and r, rho )*/
#define INTEGRATOR 1

/* 0 Heun like normal scheme
   1 Heun like in Elastics Paper (no predictor-correction for r)*/
#define HEUN_LIKE_PAPER 0

/* for simulating Solids set to 1, if 0 assuming gas*/
#define SOLIDS 0

#endif // SPH_PARAMETER_H