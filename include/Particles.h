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

#if SOLIDS
        double mu; //shear modulus
        double rhoRel; // density in relaxed state

        double *stress; // stress tensor - for each particle DIM*DIM, for all particles DIM*DIM*Number of particles

        // entries of deviatoric stress tensor - only 5 needed
        double *S12, *S13, *S23, *S11, *S22; // S21 = S21, S13 = S31, S23 = S32, S33 = -(S11+S22)

        // TO DO: add further variables needed
        
        // set parameters
        void setMu(double muSet);
        void setRho_0(double rho_0Set);
#endif


        //compute nearest neighbors in square
        void compNNSquare();
        
        // compute density of all Particles via Kernel sum
        void compDensity();
    
        /*compute change of density via continuity equation
        drho_a = rho_a sum_b m_b/rho_b * (v_a - v_b)* nabla_a W_ab
        */
        void compDrho();

        /*compute pressure of all particles, dependent of sound speed c_s
        for gases: p = c_s^2*rho 
        for solids: p = c_s^2*(rho - rho_O), rho_0 density in relaxed state*/
        void compPressure(double c_s);

        //comp acceleration of all particles
        void compAcceleration();

        // write particles to file, copied /adapted from meshless hydro
        void write2file(std::string filename);
        
};

#endif //SPH_PARTICLE_H

