#include "InitialDistribution.h" 


InitialDistribution::InitialDistribution(const std::string &file){
    HighFive::File h5file(file, HighFive::File::ReadOnly);

    // read datasets from file
    HighFive::DataSet vel = h5file.getDataSet("v");
    HighFive::DataSet mass = h5file.getDataSet("m");
    HighFive::DataSet pos = h5file.getDataSet("r");

    // read data into containers
    mass.read(m);
    pos.read(x);
    vel.read(v);

    // sanity check
    if(!(x.size() == v.size() && x.size() == m.size() && v.size() == m.size())){
        throw std::length_error("Length mismatch between velocity, position and/or mass");
    }


}

void InitialDistribution::getAllParticles(Particles &particles){
    std::vector<double>::iterator mit = m.begin();
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
        

        //std::cout << "iterator: " << (*xit)[0] << std::endl;
        //std::cout << "class variable: " << particles.y[pCounter] << std::endl;
        ++xit;
        ++vit;
        ++mit;
        ++pCounter;
        


    }

}
