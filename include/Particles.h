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
        double rhoRel; // density in relaxed state

        void setRho_0(double rho_0Set);

#if SOLIDS
        double mu; //shear modulus
        double *stress; // stress tensor - for each particle DIM*DIM, for all particles DIM*DIM*Number of particles

        // entries of deviatoric stress tensor - only 5 needed
        double *S12, *S13, *S23, *S11, *S22; // S21 = S21, S13 = S31, S23 = S32, S33 = -(S11+S22)
        // derivatives of the deviatoric stress tensor with respect to time
        double  *dS12, *dS13, *dS23, *dS11, *dS22;

        double *partialV; // partial derivatives of V

        // TO DO: add further variables needed
        // set parameters
        void setMu(double muSet);

        /*calc Stress tensor for every particle */
          
        void compStress();
        /* calc partial derivatives of v for every particle with respect to their position
           (dell v_i/ dell x_j)_a = sum_b m_b/rho_b( v_b - v_a) dell W_ab/dell x_j_a*/
        void compPartialVs();
        // calc time derivatives of deviatoric stress tensor
        void compdS();
        
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

