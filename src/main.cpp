#include <iostream>
#include "InitialDistribution.h"
#include "Integrator.h"



int N = 125; //number of particles
double sml = 2.5; // smoothing length
double c_s = 1.5; // speed of sound
double tStart = 0;
double tEnd = 1.5; // simulation end time
double deltaT = 0.01; // timestep for integration
double t = tStart; // time variable and start time




int main(int argc, char** argv){
    
    Particles particles(N, sml);
    InitialDistribution initDist("cubic_N125.h5");
    initDist.getAllParticles(particles);
    particles.compNNSquare();
    particles.compDensity();
    /* for( int i=0; i< N; i++){
        std::cout << particles.rho[i] << "\n";
        
    }
    std::cout << "\n "; */

    particles.compPressure(c_s);
    /* for( int i=0; i< N; i++){
        std::cout << particles.p[i] <<  "\n";
        
    }
    std::cout << "\n "; */
    particles.compAcceleration();
    /*     for( int i=0; i< N; i++){
        std::cout << particles.ax[i] << "  " << particles.ay[i] << "   " << particles.az[i] << "\n";
        
    }
    std::cout << "\n "; */
    particles.write2file("output/timestep0" + std::string(".h5"));
    int counter = 1;
    while ( t < tEnd){
        
        // do Timestep with integration
        std::cout << "********************+*********  time:  " << t << "   ********************" << std::endl;
        doTimestep(particles, N, sml, deltaT, c_s);
           for( int i=0; i< N; i++){
        std::cout << particles.x[i] << "  " << particles.y[i] << "   " << particles.z[i] << "\n";
        
    }
    std::cout << "\n ";
        
        if(counter%10 == 0){
     
    //write particles to file
    particles.write2file("output/timestep" +  std::to_string(counter) + std::string(".h5"));

        }
    
    //update time
    t += deltaT;
    counter++;

    }



    /* for(int i= 1; i < particles.N*particles.N; i++){
        std::cout << particles.NNsquare[i-1] << "  ";
        if(i %125 ==0){
            std::cout << "\n";
        }
    } */
    

    
   
    return 0;
}