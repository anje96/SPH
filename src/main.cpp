#include <iostream>
#include "InitialDistribution.h"


const int Dim = 3; //Dimension
int N = 125; //number of particles




int main(int argc, char** argv){
    std::cout << "Hello World" << std::endl;

    Particles particles(N);
    
    
 
    //InitialDistribution::InitialDistribution(../cubic/cubic_N125.h5);
    HighFive::File h5file("cubic_N125.h5", HighFive::File::ReadOnly);
    HighFive::DataSet vel = h5file.getDataSet("v");
    std::vector<double> v {};
    vel.read(v);
    
    
    std::vector<double>::iterator mit = v.begin();

    int pCounter =0;
    while(mit != v.end()){
        particles.m[pCounter] = *mit;
        ++mit;
        ++pCounter;
        std::cout << particles.m[pCounter] << std::endl;
    }


    

    
    

    
    return 0;
}





// int main(int argc, char** argv) {
//     std::cout << "Have " << argc << " arguments:" << std::endl;
//     for (int i = 0; i < argc; ++i) {
//         std::cout << argv[i] << std::endl;
//     }
// }