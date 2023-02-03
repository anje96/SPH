#include <iostream>
#include "InitialDistribution.h"
#include "Integrator.h"
#include "Logger.h"



int N = 125; //number of particles
double sml = 2.5; // smoothing length
double c_s = 3; // speed of sound
double tStart = 0;
double tEnd = 3; // simulation end time
double deltaT = 0.02; // timestep for integration
double t = tStart; // time variable and start time
int dumpInterval = 15; //  interval in which particles are written to file


structlog LOGCFG = {};

int main(int argc, char** argv){
    // initialize Logger
    LOGCFG.headers = true;
    LOGCFG.level = INFO;
    
   
    Logger(INFO) << " >   Reading inital distribution ... ";
    InitialDistribution initDist("cubic_N125.h5");
     Logger(INFO) << " >   Creating particles ... ";
    Particles particles(initDist.getNumberOfParticles(), sml);
    Logger(INFO) << " >   writing data to Particles...";
    initDist.getAllParticles(particles);
    Logger(INFO) << " >   calc NN in cube... ";
    particles.compNNSquare();
    Logger(INFO) << " >   calc density, pressure and acc.. ";
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
    Logger(INFO) << " >   write initial distribution to file ...";
    particles.write2file("output/timestep0" + std::string(".h5"));
    int counter = 1;
    Logger(INFO) << " >   start integration...";

    while ( t < tEnd){
        
        Logger(INFO) << " >   time: " << t;
        doTimestepHeun(particles, N, sml, deltaT, c_s);
           /* for( int i=0; i< N; i++){
        std::cout << particles.x[i] << "  " << particles.y[i] << "   " << particles.z[i] << "\n";
        
    }
    std::cout << "\n "; */
        
        if(counter%dumpInterval == 0){
     
    //write particles to file
    Logger(INFO) << " >   write particles to file...";
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