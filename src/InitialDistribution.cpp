//
// adapted from https://github.com/jammartin/meshlessHydro/blob/main/demonstrator/src/InitialDistribution.cpp
// by Anne Vera Jeschke December 2022
//
#include "InitialDistribution.h" 


InitialDistribution::InitialDistribution(const std::string &file){
    HighFive::File h5file(file, HighFive::File::ReadOnly);

    // read datasets from file
    HighFive::DataSet vel = h5file.getDataSet("v");
    HighFive::DataSet mass = h5file.getDataSet("m");
    HighFive::DataSet pos = h5file.getDataSet("x");
    
#if CALC_DENSITY == 1
    //if(h5file.exist("rho")){
        HighFive::DataSet rho = h5file.getDataSet("rho");
        
    /*}
    else{
        throw "no inital density distribution provided";
    }*/
    
#endif

    // read data into containers
    mass.read(m);
    pos.read(x);
    vel.read(v);
#if CALC_DENSITY == 1
    rho.read(density);
#endif

    // sanity check
    if(!(x.size() == v.size() && x.size() == m.size() && v.size() == m.size() 
#if CALC_DENSITY == 1
    && v.size() == density.size()
#endif
    )){
        throw std::length_error("Length mismatch between velocity, position and/or mass");
    }
    else{
        numberOfParticles = x.size();
    }


}

void InitialDistribution::getAllParticles(Particles &particles){
    std::vector<double>::iterator mit = m.begin();
#if CALC_DENSITY
    std::vector<double>::iterator rhoit = density.begin();
#endif
    std::vector<std::vector<double>>::iterator xit = x.begin();
    std::vector<std::vector<double>>::iterator vit = v.begin();
    int pCounter = 0;

    while (xit != x.end()){
        particles.m[pCounter] = *mit;
        particles.x[pCounter] = (*xit)[0];
        particles.vx[pCounter] = (*vit)[0];
        particles.y[pCounter] = (*xit)[1];
        particles.vy[pCounter] = (*vit)[1];

        particles.z[pCounter] = (*xit)[2];
        particles.vz[pCounter] = (*vit)[2];
        
#if CALC_DENSITY == 1
        particles.rho[pCounter] = *rhoit;
        rhoit++;
#endif
        //std::cout << "iterator: " << (*xit)[0] << std::endl;
        //std::cout << "class variable: " << particles.y[pCounter] << std::endl;
        ++xit;
        ++vit;
        ++mit;
        ++pCounter;
        


    }
    //printf("particle Counter in Initial Distribution: %d ", pCounter);

}
