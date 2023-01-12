#ifndef SPH_INTEGRATOR_H
#define SPH_INTEGRATOR_H

#include <iostream>
#include "parameter.h"
#include "Particles.h"

//integrate with Leap Frog Integrator
void doTimestep(Particles &particles, int numParticles, double smoothingLength, double deltaT, double c_s);


#endif // SPH_INTEGRATOR_H