//
// created by Anne Vera Jeschke December 2022
//
#include "Integrator.h"

void doTimestep(Particles &particles, double smoothingLength, double deltaT, double c_s, double rho_0, int maxNearestNeighbors){
        int numParticles = particles.N;
     // second set of particles to calc state variables and r,v,a for time t + deltaT/2
    Particles particles2(numParticles, smoothingLength, c_s, maxNearestNeighbors);

     for(int pCounter = 0; pCounter < numParticles ; pCounter++){

        //r(t + deltaT/2) = r(t)+ deltaT/2 *v(t)
        particles2.x[pCounter] = particles.x[pCounter] + deltaT*particles.vx[pCounter]/2;
        particles2.y[pCounter] = particles.y[pCounter] + deltaT*particles.vy[pCounter]/2;;
        particles2.z[pCounter] = particles.z[pCounter] + deltaT*particles.vz[pCounter]/2;;

        //v(t + deltaT/2) = v(t)+ deltaT/2 *a(t)
        particles2.vx[pCounter] = particles.vx[pCounter] + deltaT * particles.ax[pCounter]/2;
        particles2.vy[pCounter] = particles.vy[pCounter] + deltaT * particles.ay[pCounter]/2;
        particles2.vz[pCounter] = particles.vz[pCounter] + deltaT * particles.az[pCounter]/2;

        //copy masses accordingly
        particles2.m[pCounter] = particles.m[pCounter];
     }
    particles2.setRho_0(rho_0);
    // update state variable and acceleration of second set of particles
    particles2.compNNSquare();
    particles2.compDensity();

    particles2.compPressure();

    particles2.compAcceleration();


    // update main set of particles

     for(int pCounter = 0; pCounter < numParticles ; pCounter++){

        // v(t+ deltaT) = v(t) + deltaT* a(t+ deltaT/2)
        particles.vx[pCounter] = particles.vx[pCounter] + deltaT * particles2.ax[pCounter];
        particles.vy[pCounter] = particles.vy[pCounter] + deltaT * particles2.ay[pCounter];
        particles.vz[pCounter] = particles.vz[pCounter] + deltaT * particles2.az[pCounter];

        // r(t+ deltaT) = r(t+ deltaT/2) + deltaT/2* v(t+ deltaT)
        particles.x[pCounter] = particles2.x[pCounter] + deltaT*particles.vx[pCounter]/2;
        particles.y[pCounter] = particles2.y[pCounter] + deltaT*particles.vy[pCounter]/2;;
        particles.z[pCounter] = particles2.z[pCounter] + deltaT*particles.vz[pCounter]/2;;

    }

    particles.compNNSquare();
    particles.compDensity();
    particles.compPressure();
    particles.compAcceleration();
    
}


void doTimestepHeun(Particles &particles, double smoothingLength, double deltaT, double c_s, double rho_0, int maxNearestNeighbors){
    int numParticles = particles.N;
    Particles predParticles(numParticles, smoothingLength,c_s, maxNearestNeighbors);
    bool prednegativeDensity = false;

    for( int pCounter =0; pCounter < particles.N; pCounter++){
        // copy masses
        predParticles.m[pCounter] = particles.m[pCounter];
        // calc predicted velocity
        predParticles.vx[pCounter] = particles.vx[pCounter] + deltaT*particles.ax[pCounter];
        predParticles.vy[pCounter] = particles.vy[pCounter] + deltaT*particles.ay[pCounter];
        predParticles.vz[pCounter] = particles.vz[pCounter] + deltaT*particles.az[pCounter];
        

    #if HEUN_LIKE_PAPER == 1
        // calc positions (not corrected), see elastics paper
        predParticles.x[pCounter] = particles.x[pCounter] + deltaT*particles.vx[pCounter]+ 0.5*deltaT*deltaT*particles.ax[pCounter];
        predParticles.y[pCounter] = particles.y[pCounter] + deltaT*particles.vy[pCounter]+ 0.5*deltaT*deltaT*particles.ay[pCounter];
        predParticles.z[pCounter] = particles.z[pCounter] + deltaT*particles.vz[pCounter]+ 0.5*deltaT*deltaT*particles.az[pCounter];
    #else
        // calc predicted position 
        predParticles.x[pCounter] = particles.x[pCounter] + deltaT*particles.vx[pCounter];
        predParticles.y[pCounter] = particles.y[pCounter] + deltaT*particles.vy[pCounter];
        predParticles.z[pCounter] = particles.z[pCounter] + deltaT*particles.vz[pCounter];
    #endif

    predParticles.setRho_0(rho_0);
        // calc predicted density via continuity equation
#if SOLIDS
    predParticles.setMu(particles.mu);
#endif

#if CALC_DENSITY == 1
        predParticles.rho[pCounter] = particles.rho[pCounter] + deltaT* particles.drho[pCounter];
#if DEBUG_LEVEL == (1 || 2)
        
        if(predParticles.rho[pCounter]< 0){
            prednegativeDensity = true;
            #if DEBUG_LEVEL == 2
            Logger(WARN) << " Density is negative... pred particle:  " << pCounter; 
            #endif
        }
             
        
#endif  
#endif

    #if SOLIDS
    // calc predicted deviatoric stress tensor
    predParticles.S11[pCounter] = particles.S11[pCounter] + deltaT* particles.dS11[pCounter];
    predParticles.S12[pCounter] = particles.S12[pCounter] + deltaT* particles.dS12[pCounter];
    predParticles.S13[pCounter] = particles.S13[pCounter] + deltaT* particles.dS13[pCounter];
    predParticles.S22[pCounter] = particles.S22[pCounter] + deltaT* particles.dS22[pCounter];
    predParticles.S23[pCounter] = particles.S23[pCounter] + deltaT* particles.dS23[pCounter];
    #endif
   
     }
     if(prednegativeDensity){
        Logger(WARN) << " Density somewhere negative predicted Particles...";
    }

    // calc nearest neighbors
    predParticles.compNNSquare();

    // calc density via kernel sum
#if CALC_DENSITY == 0
    predParticles.compDensity();
#endif
#if CALC_DENSITY == 1

        predParticles.compDrho();
#endif

     //calc pressure
    predParticles.compPressure();

#if SOLIDS
    predParticles.compStress();
    predParticles.compPartialVs();
    predParticles.compdS();
#endif


     // calc acceleration
    predParticles.compAcceleration();

    //predParticles.write2file("output/pred" +  std::string(".h5"));



    // update particles
    bool negativeDensity = false;
    bool rIsnan = false;
    bool vIsnan = false;
    bool densityIsnan = false;
    for( int pCounter = 0; pCounter < particles.N; pCounter++){

#if HEUN_LIKE_PAPER == 1
        particles.x[pCounter] = predParticles.x[pCounter];
        particles.y[pCounter] = predParticles.y[pCounter];
        particles.z[pCounter] = predParticles.z[pCounter];
#else
        particles.x[pCounter] = predParticles.x[pCounter]+ 0.5*deltaT*(predParticles.vx[pCounter]-particles.vx[pCounter]);
        particles.y[pCounter] = predParticles.y[pCounter]+ 0.5*deltaT*(predParticles.vy[pCounter]-particles.vy[pCounter]);
        particles.z[pCounter] = predParticles.z[pCounter]+ 0.5*deltaT*(predParticles.vz[pCounter]-particles.vz[pCounter]);

#endif

        particles.vx[pCounter] = predParticles.vx[pCounter] + 0.5 *deltaT *(predParticles.ax[pCounter]-particles.ax[pCounter]);
        particles.vy[pCounter] = predParticles.vy[pCounter] + 0.5 *deltaT *(predParticles.ay[pCounter]-particles.ay[pCounter]);
        particles.vz[pCounter] = predParticles.vz[pCounter] + 0.5 *deltaT *(predParticles.az[pCounter]-particles.az[pCounter]);

#if CALC_DENSITY == 1
        particles.rho[pCounter] = predParticles.rho[pCounter] + 0.5* deltaT*(predParticles.drho[pCounter]-particles.drho[pCounter]);
        
#if DEBUG_LEVEL == (1 || 2)
        
        if(particles.rho[pCounter]< 0){
            negativeDensity = true;
            #if DEBUG_LEVEL == 2
            Logger(WARN) << " Density is negative... particle:  " << pCounter; 
            #endif
        }
        if(std::isnan(particles.rho[pCounter])){
            densityIsnan = true;
        }
        if(std::isnan(particles.x[pCounter]) || std::isnan(particles.y[pCounter])||std::isnan(particles.y[pCounter])){
            rIsnan = true;
        }
         if(std::isnan(particles.vx[pCounter]) || std::isnan(particles.vy[pCounter])||std::isnan(particles.vy[pCounter])){
            vIsnan = true;
        }
        
        
#endif
#endif

#if SOLIDS
        particles.S11[pCounter] = predParticles.S11[pCounter] + 0.5* deltaT*(predParticles.dS11[pCounter]- particles.dS11[pCounter]);
        particles.S12[pCounter] = predParticles.S12[pCounter] + 0.5* deltaT*(predParticles.dS12[pCounter]- particles.dS12[pCounter]);
        particles.S13[pCounter] = predParticles.S13[pCounter] + 0.5* deltaT*(predParticles.dS13[pCounter]- particles.dS13[pCounter]);
        particles.S22[pCounter] = predParticles.S22[pCounter] + 0.5* deltaT*(predParticles.dS22[pCounter]- particles.dS22[pCounter]);
        particles.S23[pCounter] = predParticles.S23[pCounter] + 0.5* deltaT*(predParticles.dS23[pCounter]- particles.dS23[pCounter]);
#endif
 
    }
    if(negativeDensity){
        Logger(WARN) << " Density somewhere negative...";
    }
    if(densityIsnan){
        Logger(WARN) << " rho nan somewhere...";
    }
    if(rIsnan){
        Logger(WARN) << " r nan somewhere...";
    }
     if(vIsnan){
        Logger(WARN) << " v nan somewhere...";
    }
    
    particles.compNNSquare();
    
#if CALC_DENSITY == 0
    particles.compDensity();
#else // calc drho for integration

    particles.compDrho();
#endif

    particles.compPressure();

#if SOLIDS
    particles.compStress();
    particles.compPartialVs();
    particles.compdS();
#endif

    particles.compAcceleration();    

}







