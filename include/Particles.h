//
// Created by Anne Vera Jeschke December 2022
//
#ifndef SPH_PARTICLES_H
#define SPH_PARTICLES_H

#include <iostream>
#include <cmath>
#include <limits>
#include "Kernel.h"
#include <highfive/H5File.hpp>
#include "Logger.h"





class Particles{
    public: 
        Particles(int NumParticles, double smoothingLength, double speedOfSound, int maxNearestNeighbors);

        ~Particles();

        int N, maxNN; // maxNN is the max number of nearest neighbors
        double sml;
        double *m, *x, *y, *z, *vx, *vy, *vz, *rho, *drho, *p, *ax, *ay, *az;
        int *iCounter; //counts interactions of the particles
        int *NNsquare;
        double rhoRel; // density in relaxed state
        double c_s; // speed of sound
        

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
        int indexSigma(int pCounter, int i , int j);
        int indexMatrix(int i, int j);

        /*calc Stress tensor for every particle */
          
        void compStress();
        /* calc partial derivatives of v for every particle with respect to their position
           (dell v_i/ dell x_j)_a = sum_b m_b/rho_b( v_b - v_a) dell W_ab/dell x_j_a*/
        void compPartialVs();
        // calc time derivatives of deviatoric stress tensor
        void compdS();
        
#endif


#if ARTIFICIAL_VISCOSITY 
        /* free parameters of artificial viscosity*/
        double alpha = 1; 
        double beta = 2*alpha;

        // no divergence for artificial viscosity
        double epsilon = 0.01;

        // calc artificial viscosity term for particle pCounter and neighbor
        double compArtificialVisc( int pCounter, int neighbor);
#endif

#if ARTIFICIAL_STRESS
        //epsilon of artificial stress
        double epsilonAS = 0.3;
        // particle spacing
        double deltaP = 0.1;
        // for artificial stress term R_ab * f^n
        int n = 6;
        void compRij(int pCounter, int neighbor, double *R_ab);
        double compFn(int pCounter, int neighbor);
#endif

#if COURANT_CONDITION
    double compTimestep(double timestep);
    double CFL = 0.4; //CFL-Number
    double mu_max;
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
        void compPressure();

        //comp acceleration of all particles
        void compAcceleration();

        // write particles to file, copied /adapted from meshless hydro
        void write2file(std::string filename);
        
};

#endif //SPH_PARTICLE_H

