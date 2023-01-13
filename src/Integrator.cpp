#include "Integrator.h"

void doTimestep(Particles &particles, int numParticles, double smoothingLength, double deltaT, double c_s){

     // second set of particles to calc state variables and r,v,a for time t + deltaT/2
    Particles particles2(numParticles, smoothingLength);

     for(int pCounter = 0; pCounter < numParticles ; pCounter++){

        //r(t + deltaT/2) = r(t)+ deltaT/2 *v(t)
        particles2.x[pCounter] = particles.x[pCounter] + deltaT*particles.vx[pCounter]/2;
        particles2.y[pCounter] = particles.y[pCounter] + deltaT*particles.vy[pCounter]/2;;
        particles2.z[pCounter] = particles.z[pCounter] + deltaT*particles.vz[pCounter]/2;;

         /*    for( int i=0; i< numParticles; i++){
        std::cout << particles2.x[i] << "  " << particles2.y[i] << "   " << particles2.z[i] << "\n";
        
        }
        std::cout << "\n ";
 */
        //v(t + deltaT/2) = v(t)+ deltaT/2 *a(t)
        particles2.vx[pCounter] = particles.vx[pCounter] + deltaT * particles.ax[pCounter]/2;
        particles2.vy[pCounter] = particles.vy[pCounter] + deltaT * particles.ay[pCounter]/2;
        particles2.vz[pCounter] = particles.vz[pCounter] + deltaT * particles.az[pCounter]/2;

        /*  for( int i=0; i< numParticles; i++){
        std::cout << particles2.vx[i] << "  " << particles2.vy[i] << "   " << particles2.vz[i] << "\n";
        
        }
        std::cout << "\n "; */

        //copy masses accordingly
        particles2.m[pCounter] = particles.m[pCounter];
     }

    // update state variable and acceleration of second set of particles
    particles2.compNNSquare();
    particles2.compDensity();
   /*  std::cout << "Density: " << std::endl;
    for( int i=0; i< numParticles; i++){
        std::cout << particles2.rho[i] << "\n";
        
    }
    std::cout << "\n "; */

    particles2.compPressure(c_s);

    /* std::cout << "Pressure: " << std::endl;
    for( int i=0; i< numParticles; i++){
        std::cout << particles2.p[i] <<  "\n";
        
    }
    std::cout << "\n "; */
    particles2.compAcceleration();

    /* std::cout << "Acceleration: " << std::endl;
        for( int i=0; i< numParticles; i++){
        std::cout << particles2.ax[i] << "  " << particles2.ay[i] << "   " << particles2.az[i] << "\n";
        
    }
    std::cout << "\n "; */

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
    particles.compPressure(c_s);
    particles.compAcceleration();
    
}


