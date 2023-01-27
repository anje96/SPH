#ifndef SPH_PARTICLES_H
#define SPH_PARTICLES_H

#include <iostream>
#include "Kernel.h"
#include <highfive/H5File.hpp>
#include "Logger.h"




class Particles{
    public: 
        Particles(int NumParticles, double smoothingLength);
        ~Particles();

        int N;
        double sml;
        double *m, *x, *y, *z, *vx, *vy, *vz, *rho, *drho, *p, *ax, *ay, *az;
        int *NNsquare;

        //compute nearest neighbors in square
        void compNNSquare();
        
        // compute density of all Particles via Kernel sum
        void compDensity();
    
        // compute change of density via continuity equation
        void ChangeOfDensity();

        //compute pressure of all particles, dependent of sound speed
        void compPressure(double c_s);

        //comp acceleration of all particles
        void compAcceleration();

        // write particles to file, copied /adapted from meshless hydro
        void write2file(std::string filename);
        
};

#endif //SPH_PARTICLE_H

