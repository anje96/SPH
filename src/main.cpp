#include <iostream>
#include "InitialDistribution.h"
#include "Integrator.h"
#include "Logger.h"
#include <cxxopts.hpp>
#include "ConfigParser.h"


structlog LOGCFG = {};

int main(int argc, char** argv){

    cxxopts::Options cmdLineOptions { "sph",
                                      "testcode for isothermal gas and solids" };
    cmdLineOptions.add_options()
            ("c,config", "Path to config file", cxxopts::value<std::string>()->default_value("config.info"))
          /*   ("v,verbose", "More printouts for debugging")
            ("s,silent", "Suppress normal printouts")
            ("h,help", "Show this help") */;

    auto cmdLineOpts = cmdLineOptions.parse(argc, argv);

  /*   if (cmdLineOpts.count("help")) {
        std::cout << cmdLineOptions.help() << std::endl;
        exit(0);
    } */

    ConfigParser confP { cmdLineOpts["config"].as<std::string>() };

    // initialize Logger
    LOGCFG.headers = true;
    LOGCFG.level = INFO;

    Logger(INFO) << "Reading configuration ... ";
    ConfigParser::Configuration config; // initialize configuration

    config.initFile = confP.getVal<std::string>("initFile");  
    Logger(INFO) << "   >   Initial distribution: " << config.initFile;
    config.outDir = confP.getVal<std::string>("outDir");
    Logger(INFO) << "   >   Output directory: " << config.outDir;
    config.timeStep = confP.getVal<double>("timeStep");
    Logger(INFO) << "   >   Time step: " << config.timeStep;
    config.timeEnd = confP.getVal<double>("timeEnd");
    Logger(INFO) << "   >   End of simulation: " << config.timeEnd;
    config.h5DumpInterval = confP.getVal<int>("h5DumpInterval");
    Logger(INFO) << "   >   Dump data to h5 file every " << config.h5DumpInterval << " steps";
    config.speedOfSound = confP.getVal<double>("speedOfSound");
    Logger(INFO) << "   >   Speed of Sound: " << config.speedOfSound;
    config.smoothingLength = confP.getVal<double>("smoothingLength");
    Logger(INFO) << "   >   smoothing Length: " << config.smoothingLength;
    config.restDensity = confP.getVal<double>("restDensity");
    Logger(INFO) << "   >   Density of particles in relaxed State, should be zero for gas: " << config.restDensity;
    config.maxNN = confP.getVal<int>("maxNN");
    Logger(INFO) << "   >   maximum number of nearest neighbors: " << config.maxNN;

#if SOLIDS
    config.shearModulus = confP.getVal<double>("shearModulus");
    Logger(INFO) << "   >   Shear modulus: " << config.shearModulus;
#endif

    Logger(INFO) << "   >   Reading initial distribution ... ";
    InitialDistribution initDist(config.initFile);
    Logger(INFO) << "   >   Creating particles ... ";
    Particles particles(initDist.getNumberOfParticles(), config.smoothingLength, config.speedOfSound, config.maxNN);
    Logger(INFO) << "   >   number of particles: "<< initDist.getNumberOfParticles();
    Logger(INFO) << "   >   writing data to Particles...";
    initDist.getAllParticles(particles);

#if SOLIDS
    Logger(INFO) << "   >   setting shear module...";
    particles.setMu(config.shearModulus);
    // TO DO: add further calls when necessary
#endif
    Logger(INFO) << "   >   setting rest density:...";
    particles.setRho_0(config.restDensity); // set rho_0 to zero for gas



    Logger(INFO) << "   >   calc NN in cube... ";
    particles.compNNSquare();


    

#if CALC_DENSITY == 0 // Kernel sum
    Logger(INFO) << "   >   calc density... ";
    particles.compDensity();
#else // calc drho for later integration of rho
    Logger(INFO) << "   >   calc drho... ";
    particles.compDrho();
    //particles.drho[0] = 1;
#endif
   
    Logger(INFO) << "   >   calc pressure...";
    particles.compPressure();
    //particles.p[0] = 2;

#if SOLIDS
    Logger(INFO) << "   >   calc stress...";
    particles.compStress();


    Logger(INFO) << "   >   calc partial derivatives of v with respect to x...";
    particles.compPartialVs();

    Logger(INFO) << "   >   calc dS/dt...";
    particles.compdS();
#endif

    Logger(INFO) << "   >   calc acc...";
    particles.compAcceleration();
   
    Logger(INFO) << "   >   write initial distribution to file ...";
    particles.write2file(config.outDir +"/timestep0000" + std::string(".h5"));
    int counter = 1;
    Logger(INFO) << "   >   start integration...";

    double t = 0;
    int writeCounter = 0;

    while ( t < config.timeEnd){
        
        Logger(INFO) << "   >   time: " << t;
        doTimestepHeun(particles, config.smoothingLength, config.timeStep, config.speedOfSound, config.restDensity, config.maxNN);        
        if(counter%config.h5DumpInterval == 0){
     
    //write particles to file
    Logger(INFO) << "   >   write particles to file... timestep: " << counter;
    std::string number = std::to_string(writeCounter);
    number.insert(number.begin(), 4 - number.length(), '0');
    
    particles.write2file(config.outDir+"/timestep" + number + std::string(".h5"));
    writeCounter++;
        }
    
    //update time
    t += config.timeStep;
    counter++;

    }
   
   
    return 0;
}