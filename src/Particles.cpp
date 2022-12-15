#include "../include/Particles.h" 

//int NumParticles = 125;
Particles::Particles( int NumParticles){
    x = new double[NumParticles];
    y = new double[NumParticles];
    z = new double[NumParticles];
    vx = new double[NumParticles];
    vy = new double[NumParticles];
    vz = new double[NumParticles];
    m = new double[NumParticles];
    rho = new double[NumParticles];
    p = new double[NumParticles];
    ax = new double[NumParticles];
    ay = new double[NumParticles];
    az = new double[NumParticles];


};

Particles::~Particles(){
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] vx;
    delete[] vy;
    delete[] vz;
    delete[] m;
    delete[] rho;
    delete[] p;
    delete[] ax;
    delete[] ay;
    delete[] az;

}
