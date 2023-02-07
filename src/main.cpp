#include <iostream>
#include "InitialDistribution.h"
#include "Integrator.h"
#include "Logger.h"
#include <cxxopts.hpp>
#include "ConfigParser.h"


structlog LOGCFG = {};

int main(int argc, char** argv){

    cxxopts::Options cmdLineOptions { "mlh",
                                      "Demonstrator for the meshless hydrodynamic simulation methods MFV and MFM" };
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

#if SOLIDS
    config.restDensity = confP.getVal<double>("restDensity");
    Logger(INFO) << "   >   Density of particles in relaxed State: " << config.restDensity;
    config.shearModulus = confP.getVal<double>("shearModulus");
    Logger(INFO) << "   >   Shear modulus: " << config.shearModulus;
#endif

       Logger(INFO) << "   >   Reading inital distribution ... ";
    InitialDistribution initDist(config.initFile);
     Logger(INFO) << "   >   Creating particles ... ";
    Particles particles(initDist.getNumberOfParticles(), config.smoothingLength);
    Logger(INFO) << "   >   writing data to Particles...";
    initDist.getAllParticles(particles);

#if SOLIDS
    Logger(INFO) << "   >   setting shear module, rest density...";
    particles.setMu(config.shearModulus);
    particles.setRho_0(config.restDensity);
    // TO DO: add further calls when necessary
#endif

    Logger(INFO) << "   >   calc NN in cube... ";
    particles.compNNSquare();

    Logger(INFO) << "   >   calc density... ";

#if CALC_DENSITY == 0 // Kernel sum
    particles.compDensity();
#else // calc drho for later integration of rho
    particles.compDrho();
#endif

   
    Logger(INFO) << "   >   calc pressure and acc.. ";
    particles.compPressure(config.speedOfSound);
    particles.compAcceleration();
   
    Logger(INFO) << "   >   write initial distribution to file ...";
    particles.write2file(config.outDir +"/timestep0" + std::string(".h5"));
    int counter = 1;
    Logger(INFO) << "   >   start integration...";

    double t = 0;

    while ( t < config.timeEnd){
        
        Logger(INFO) << "   >   time: " << t;
        doTimestepHeun(particles, config.smoothingLength, config.timeStep, config.smoothingLength);        
        if(counter%config.h5DumpInterval == 0){
     
    //write particles to file
    Logger(INFO) << "   >   write particles to file...";
    particles.write2file(config.outDir+"/timestep" +  std::to_string(counter) + std::string(".h5"));

        }
    
    //update time
    t += config.timeStep;
    counter++;

    }
   
   
    return 0;
}