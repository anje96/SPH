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
    drho = new double[NumParticles];
    p = new double[NumParticles];
    ax = new double[NumParticles];
    ay = new double[NumParticles];
    az = new double[NumParticles]; 
    
    N = NumParticles;
    sml = smoothingLength;

#if SOLIDS

    /* access for \sigma_ij of particle p (particle index 0 to N-1) via stress[p*DIM*DIM + i*DIM + j]
    */
    stress = new double[NumParticles*DIM*DIM]; 
    S11 = new double[NumParticles];
    S12 = new double[NumParticles];
    S13 = new double[NumParticles];
    S22 = new double[NumParticles];
    S23 = new double[NumParticles];
    dS11 = new double[NumParticles];
    dS12 = new double[NumParticles];
    dS13 = new double[NumParticles];
    dS22 = new double[NumParticles];
    dS23 = new double[NumParticles];
    /* access for dV_i/dx_j (i row, j column) of particle p (particle index 0 to N-1) via stress[p*DIM*DIM + i*DIM + j]
    */
    partialV = new double[NumParticles*DIM*DIM];

#endif

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
    delete[] drho;
    delete[] p;
    delete[] ax;
    delete[] ay;
    delete[] az; 

#if SOLIDS
    delete[] stress;
    delete[] S11;
    delete[] S12;
    delete[] S13;
    delete[] S22;
    delete[] S23;
    delete[] dS11;
    delete[] dS12;
    delete[] dS13;
    delete[] dS22;
    delete[] dS23;
    delete[] partialV;
#endif

};

// setting functions
#if SOLIDS
void Particles::setMu(double muSet){
    mu = muSet;
};
#endif

void Particles::setRho_0(double rho_0Set){
   rhoRel = rho_0Set;
};

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

        // self density of particle
        double rhoTemp = m[pCounter]*cubicSpline(x[pCounter], y[pCounter], z[pCounter], x[pCounter], y[pCounter], z[pCounter], sml);

        neighbor = NNsquare[N*pCounter+numberOfNN];
        while( neighbor != -1){
            
            rhoTemp += m[neighbor]*cubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml);
            numberOfNN++;
            neighbor = NNsquare[N*pCounter+numberOfNN];
        }
        if(rhoTemp < 0){
            Logger(WARN) << "Negative density, particle number: " << pCounter; 
        }
        rho[pCounter] = rhoTemp;        
               
    }
 
};

void Particles::compDrho(){
    for( int pCounter =0; pCounter < N; pCounter++){
        int numberOfNN = 0;
        int neighbor = 0;

        // set change to zero
        double drhoTemp = 0.0;
        double* gradKernel;
        double deltaVx;
        double deltaVy;
        double deltaVz;

        neighbor = NNsquare[N*pCounter+numberOfNN];
        while( neighbor != -1){

            gradKernel = gradCubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml);
            deltaVx = vx[pCounter]- vx[neighbor];
            deltaVy = vy[pCounter]- vy[neighbor];
            deltaVz = vz[pCounter]- vz[neighbor];
            
            drhoTemp += m[neighbor]*(deltaVx*gradKernel[0] + deltaVy*gradKernel[1] + deltaVz*gradKernel[2])/rho[neighbor];
            numberOfNN++;
            neighbor = NNsquare[N*pCounter+numberOfNN];
        }
        
        drho[pCounter] = rho[pCounter]*drhoTemp;        
               
    }
    
};

// Calc pressure over EOS, here isothermal, p dependent on sound velocity c_s
void Particles::compPressure(double c_s){
    for( int pCounter = 0; pCounter < N; pCounter++){
        p[pCounter] = c_s*c_s*(rho[pCounter] - rhoRel);
    }
    
};

#if SOLIDS
    // Calc stress tensor for every particle
    void Particles::compStress(){
        for( int pCounter = 0;  pCounter < N; pCounter++){
            double S_ij[DIM*DIM] = {S11[pCounter], S12[pCounter], S13[pCounter], S12[pCounter], S22[pCounter], S23[pCounter]
            , S13[pCounter], S23[pCounter], -(S11[pCounter]+S22[pCounter]) };
            for( int i = 0; i < DIM; i++){
                for(int j = 0; j < DIM; j++){
                    if(i == j){
                        stress[pCounter*DIM*DIM + i*DIM+ j] = -p[pCounter] + S_ij[i*DIM+j];
                    }
                    else{
                        stress[pCounter*DIM*DIM + i*DIM+ j] =  S_ij[i*DIM+j];
                    }
                   
                    
                }
            }
       }
   
    };



void Particles::compPartialVs(){
    // set every partial derivative to zero first
    for(int k = 0; k < DIM*DIM*N; k++){
        partialV[k] = 0;
    }
    // loop over particles
    for( int pCounter = 0; pCounter < N; pCounter++){
        int numberOfNN = 0;
        int neighbor = 0;
        double* gradKernel;

        double vParticle[DIM]  = {vx[pCounter], vy[pCounter], vz[pCounter]};

        neighbor = NNsquare[N*pCounter+numberOfNN];
        
        // loop over nearest neighbors of each particle to calc dvi/dxj
        while(neighbor != -1){
            gradKernel = gradCubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml);
            double vNeighbor[DIM] = {vx[neighbor], vy[neighbor], vz[neighbor]};

            for( int i = 0; i < DIM; i++){
                for( int j = 0; j < DIM; j++){
                partialV[DIM*DIM*pCounter + i*DIM + j] += m[neighbor]/rho[neighbor] *( vNeighbor[i]-vParticle[i]) * gradKernel[j];                    
                }    
            }

            // switch to  next neighbor
            numberOfNN++;
            neighbor = NNsquare[N*pCounter+numberOfNN];

        } 
        // go to next particle

    }
};


void Particles::compdS(){
    for(int pCounter = 0; pCounter < N; pCounter++){
        
        double S_ij[DIM*DIM] = {S11[pCounter], S12[pCounter], S13[pCounter], S12[pCounter], S22[pCounter], S23[pCounter]
            , S13[pCounter], S23[pCounter], -(S11[pCounter]+S22[pCounter]) }; // deviatoric stress tensor

        double epsilon[DIM*DIM]; // strain rate tensor
        double omega[DIM*DIM]; // rotation rate tensor

        for(int i = 0; i < DIM; i++){
            for(int j = 0; j < DIM; j++){
                epsilon[i*DIM+ j] = 0.5*(partialV[pCounter*DIM*DIM + i*DIM + j]+partialV[pCounter*DIM*DIM + j*DIM + i] );
                omega[i*DIM+ j] = 0.5*(partialV[pCounter*DIM*DIM + i*DIM + j]- partialV[pCounter*DIM*DIM + j*DIM + i] );
            }
        }

        dS11[pCounter] = 2*mu* (epsilon[0*DIM+0]- (1/3)* epsilon[0*DIM+0]);
        dS12[pCounter] = 2*mu* epsilon[0*DIM+1];
        dS22[pCounter] = 2*mu* (epsilon[1*DIM+1]- (1/3)* epsilon[1*DIM+1]);
        dS13[pCounter] = 2*mu* epsilon[0*DIM+2];
        dS23[pCounter] = 2*mu* epsilon[1*DIM+2];

        for(int k = 0; k < DIM; k++){
            
            dS11[pCounter] += S_ij[0*DIM+k]*omega[0*DIM+k] + omega[0*DIM+k]* S_ij[k*DIM+0];
            dS12[pCounter] += S_ij[0*DIM+k]*omega[1*DIM+k] + omega[0*DIM+k]* S_ij[k*DIM+1];
            dS22[pCounter] += S_ij[1*DIM+k]*omega[1*DIM+k] + omega[1*DIM+k]* S_ij[k*DIM+1];
            dS13[pCounter] += S_ij[0*DIM+k]*omega[2*DIM+k] + omega[0*DIM+k]* S_ij[k*DIM+2];
            dS23[pCounter] += S_ij[1*DIM+k]*omega[2*DIM+k] + omega[1*DIM+k]* S_ij[k*DIM+2];
        }


    }

};

#endif

// Calc acceleration through differential equation, iterate over nearest neighbors
// TO DO: self-acceleration? derivative for kernel zero for r_ij = 0. --> not necessary
void Particles::compAcceleration(){
    // iterate over all particles
    for( int pCounter =0; pCounter < N; pCounter++){
        int numberOfNN = 0;
        int neighbor = 0;

        double axTemp = 0;
        double ayTemp = 0;
        double azTemp = 0;

        

        double* gradKernel;
#if SOLIDS
        double dV[DIM] = {0,0,0};
#else
        double prefactor = 0;
#endif

        neighbor = NNsquare[N*pCounter+numberOfNN];
        // sum over nearest neighbors (while- loop)
        while( neighbor != -1){
#if SOLIDS
            gradKernel = gradCubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml);
            // calc dv_i for one neighbor with sum over stress tensor
            for( int i = 0; i < DIM; i++){
                for( int j = 0; j < DIM; j++){
                    dV[i] += (stress[pCounter*DIM*DIM + i*DIM+ j]/(rho[pCounter]*rho[pCounter]) + 
                    stress[neighbor*DIM*DIM + i*DIM+ j]/( rho[neighbor]*rho[neighbor]))*gradKernel[j];
                }
            }

            // add up components of the sum to the total sum for each nearest neighbor
            axTemp += m[neighbor]* dV[0];
            ayTemp += m[neighbor]* dV[1];
            azTemp += m[neighbor]* dV[2];


#else
            /*1st Version of SPH equation*/
            //prefactor = - m[neighbor] *(p[neighbor]+p[pCounter])/(rho[pCounter]*rho[neighbor]);

            /* 2nd Version of SPH equation*/
            prefactor = -m[neighbor]* (p[neighbor]/(rho[neighbor]*rho[neighbor]) + p[pCounter]/(rho[pCounter]*rho[pCounter]));
            
            
            gradKernel = gradCubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml);
            
            
            axTemp += prefactor*gradKernel[0];
            ayTemp += prefactor*gradKernel[1];
            azTemp += prefactor*gradKernel[2];
#endif           
            numberOfNN++;
            neighbor = NNsquare[N*pCounter+numberOfNN];
        }
        if(std::isnan(axTemp) || std::isnan(ayTemp) || std::isnan(azTemp) ) {
                Logger(WARN) <<" >  acceleration is nan, particle: " << pCounter;
            }
        // write new accelerations 
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

