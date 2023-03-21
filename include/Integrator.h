#ifndef SPH_INTEGRATOR_H
#define SPH_INTEGRATOR_H

#include <iostream>
#include "parameter.h"
#include "Particles.h"
#include "InitialDistribution.h"

//integrate with Leap Frog Integrator
void doTimestep(Particles &particles, double smoothingLength, double deltaT, double c_s, double rho_0, int maxNearestNeighbors);

// integrate with predrictor-corrector step, like in Elastics paper
void doTimestepHeun(Particles &particles,  double smoothingLength, double deltaT, double c_s, double rho_0, int maxNearestNeighbors);


#endif // SPH_INTEGRATOR_H