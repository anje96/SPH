#ifndef SPH_INITIALDISTRIBUTION_H
#define SPH_INITIALDISTRIBUTION_H

#include <highfive/H5File.hpp>

#include "Particles.h"


//copied or adapted from meshlesshydro
class InitialDistribution {
public:
    InitialDistribution(const std::string &file);
    void getAllParticles(Particles &particles);

private:
    //containers to be filled from hdf5 file
    std::vector<double> m {};
    std::vector<std::vector<double>> x {}, v{};
};

#endif // SPH_INITIALDISTRIBUTION_H
