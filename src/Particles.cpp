#include "../include/Particles.h" 
#include "Kernel.h"

//int NumParticles = 125;
Particles::Particles( int NumParticles, double smoothingLength) {
    // default initialization is zero
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
    
    N = NumParticles;
    sml = smoothingLength;

    /* to store index of nearest neighbors of each particle
        acces via index = NumParticles * X + Y
        X indicating particle index and Y stores the index of nearest neighbor particle, if -1 no nearest neighbor
    */
    NNsquare= new int[NumParticles*NumParticles];

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

// find nearest neighbor of each particle in cube with side length equal to smoothing length and write them into NNsquare
void Particles::compNNSquare(){
     // initialize nearest neighbors with -1
    for(int i =0; i < N*N; i++){
        NNsquare[i] = -1;
    }

    double xmin = 0.0;
    double xmax = 0.0;
    double ymin = 0.0;
    double ymax = 0.0;
    double zmin = 0.0;
    double zmax = 0.0;
    int NNCounter;
    

    for( int pCounter = 0; pCounter < N; pCounter++) {
        NNCounter = 0;
        xmin = x[pCounter]-sml;
        xmax = x[pCounter]+sml;
        ymin = y[pCounter]-sml;
        ymax = y[pCounter]+sml;
        zmin = z[pCounter]-sml;
        zmax = z[pCounter]+sml;

        for( int nCounter = 0; nCounter < N; nCounter++){
            if(nCounter == pCounter){
                continue;
            }
            else{
                if((xmin < x[nCounter] && x[nCounter] < xmax) && (ymin < y[nCounter] && y[nCounter] < ymax) && (zmin < z[nCounter] && z[nCounter] < zmax)){
                    NNsquare[pCounter * N + NNCounter] = nCounter;
                    NNCounter++;

                }

            }
            

        }
   }
};

// Calc density for all particles
void Particles::compDensity(){
      for( int pCounter =0; pCounter < N; pCounter++){
        int numberOfNN = 0;
        int neighbor = 0;

        double rhoTemp = 0;

        neighbor = NNsquare[N*pCounter+numberOfNN];
        while( neighbor != -1){
            
            rhoTemp += m[neighbor]*cubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml);
            numberOfNN++;
            neighbor = NNsquare[N*pCounter+numberOfNN];
        }
        rho[pCounter] = rhoTemp;        
               
    }
 
};

// Calc pressure over EOS, here isothermal, p dependent on sound velocity c_s
void Particles::compPressure(double c_s){
    for( int pCounter = 0; pCounter < N; pCounter++){
        p[pCounter] = c_s*c_s*rho[pCounter];
    }
    
};

// Calc acceleration through differential equation, iterate over nearest neighbors
void Particles::compAcceleration(){
    for( int pCounter =0; pCounter < N; pCounter++){
        int numberOfNN = 0;
        int neighbor = 0;

        double axTemp = 0;
        double ayTemp = 0;
        double azTemp = 0;

        double prefactor = 0;

        double* gradKernel;

        neighbor = NNsquare[N*pCounter+numberOfNN];
        while( neighbor != -1){

            //problematic if dividing through 0
            if(rho[pCounter]*rho[neighbor]== 0){
                prefactor = 0;
            }
            else{
                prefactor = - m[neighbor] *(p[neighbor]+p[pCounter])/(rho[pCounter]*rho[neighbor]);
            }
            
            gradKernel = gradCubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml);
            
            axTemp += prefactor*gradKernel[0];
            ayTemp += prefactor*gradKernel[1];
            azTemp += prefactor*gradKernel[2];
            
            numberOfNN++;
            neighbor = NNsquare[N*pCounter+numberOfNN];
        }
        ax[pCounter] = axTemp;
        ay[pCounter] = ayTemp;
        az[pCounter] = azTemp;       
               
    }

};

void Particles::write2file(std::string filename){
    // open output file
    HighFive::File h5File { filename, HighFive::File::ReadWrite |
                                      HighFive::File::Create |
                                      HighFive::File::Truncate };
    
    // dimensions for datasets containing vectors
    std::vector<size_t> dataSpaceDims(2);
    dataSpaceDims[0] = std::size_t(N); // number of particles
    dataSpaceDims[1] = DIM;

    // create datasets
    HighFive::DataSet rhoDataSet = h5File.createDataSet<double>("/rho", HighFive::DataSpace(N));
    HighFive::DataSet pDataSet = h5File.createDataSet<double>("/p", HighFive::DataSpace(N));
    HighFive::DataSet mDataSet = h5File.createDataSet<double>("/m", HighFive::DataSpace(N));
    //HighFive::DataSet uDataSet = h5File.createDataSet<double>("/u", HighFive::DataSpace(N));
    HighFive::DataSet posDataSet = h5File.createDataSet<double>("/r", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet velDataSet = h5File.createDataSet<double>("/v", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet accDataSet = h5File.createDataSet<double>("/a", HighFive::DataSpace(dataSpaceDims));

      // containers for particle data
    std::vector<double> rhoVec(rho,rho+N);
    std::vector<double> pVec(p,p+N);
    std::vector<double> mVec(m, m + N);

    std::vector<std::vector<double>> posVec(N);
    std::vector<std::vector<double>> velVec(N);
    std::vector<std::vector<double>> accVec(N);

     // fill containers with data
    std::vector<double> posBuf(DIM);
    std::vector<double> velBuf(DIM);
    std::vector<double> accBuf(DIM);
    for(int i=0; i<N; ++i){

        //position
        posBuf[0] = x[i];
        posBuf[1] = y[i];
#if DIM == 3
        posBuf[2] = z[i];
#endif
        posVec[i] = posBuf;

        // velocity
        velBuf[0] = vx[i];
        velBuf[1] = vy[i];
#if DIM == 3
        velBuf[2] = vz[i];
#endif
        velVec[i] = velBuf;

        // density gradient
        accBuf[0] = ax[i];
        accBuf[1] = ay[i];
#if DIM == 3
        accBuf[2] = az[i];
#endif
        accVec[i] = accBuf;
    }
    // write data
    rhoDataSet.write(rhoVec);
    pDataSet.write(pVec);
    mDataSet.write(mVec);
    
    posDataSet.write(posVec);
    velDataSet.write(velVec);
    accDataSet.write(accVec);   

};

