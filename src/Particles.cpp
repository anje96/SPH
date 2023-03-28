#include "Particles.h" 
#include "Kernel.h"

Particles::Particles( int NumParticles, double smoothingLength, double speedOfSound, int maxNearestNeighbors) {
   
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

    iCounter = new int[NumParticles]; // counts interactions of the particles which is equal to the number of nearest neighbors of the particle

    maxNN = maxNearestNeighbors;
    
    N = NumParticles;
    sml = smoothingLength;
    c_s = speedOfSound;

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

#if COURANT_CONDITION
    mu_max = 0.0;
#endif

    /* to store index of nearest neighbors of each particle
        access via index = NumParticles * X + Y
        X indicating particle index and Y stores the index of nearest neighbor particle, if -1 no nearest neighbor
    */
    NNsquare= new int[NumParticles * maxNN];

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
    delete[] NNsquare;

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

// returns corresponding indices for sigma_ij of particle pCounter
// or dv_i/dx_j for particle pCounter 
int Particles::indexSigma(int pCounter, int i , int j){
    return pCounter*DIM*DIM+i*DIM+j;
}
//returns corresponding indices for entry ij of a matrix of Dimension DIM*DIM 
int Particles::indexMatrix(int i, int j){
    return i*DIM+j;
}
#endif

void Particles::setRho_0(double rho_0Set){
   rhoRel = rho_0Set;
};

//helper functions



// find nearest neighbor of each particle in cube with side length equal to smoothing length and write them into NNsquare
void Particles::compNNSquare(){
     // initialize nearest neighbors with -1
    for(int i =0; i < maxNN*N; i++){
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
        NNCounter = 0; //nearest neighbor counter which are within smoothing length
        xmin = x[pCounter]-sml;
        xmax = x[pCounter]+sml;
        ymin = y[pCounter]-sml;
        ymax = y[pCounter]+sml;
        zmin = z[pCounter]-sml;
        zmax = z[pCounter]+sml;

        for( int nCounter = 0; nCounter < N && NNCounter < maxNN; nCounter++){
            if(nCounter == pCounter){
                continue;
            }
            else{
                if((xmin <= x[nCounter] && x[nCounter] <= xmax) && (ymin <= y[nCounter] && y[nCounter] <= ymax) && (zmin <= z[nCounter] && z[nCounter] <= zmax)){
                    NNsquare[pCounter * maxNN + NNCounter] = nCounter;
                    ++NNCounter;

                }

            }
            
        }

        iCounter[pCounter] = NNCounter;

   }

#if DEBUG_LEVEL == 2
    Logger(DEBUG) << "   >  Nearest Neighbors: ";
    for( int pCounter = 0; pCounter < N; pCounter++){
        Logger(DEBUG) << "Neighbors of particle: "<< pCounter << "    :";
        int numberOfNN = 0;
        int neighbor = 0;
        neighbor = NNsquare[maxNN * pCounter + numberOfNN];
        while( neighbor != -1 && numberOfNN < maxNN){
            Logger (DEBUG) << neighbor;
            numberOfNN++;
            neighbor = NNsquare[maxNN * pCounter + numberOfNN];
        }
        
    }
#endif
};

// Calc density for all particles
void Particles::compDensity(){
      for( int pCounter =0; pCounter < N; pCounter++){
        int numberOfNN = 0;
        int neighbor = 0;

        // self density of particle
        double rhoTemp = m[pCounter]*cubicSpline(x[pCounter], y[pCounter], z[pCounter], x[pCounter], y[pCounter], z[pCounter], sml);

        neighbor = NNsquare[maxNN * pCounter + numberOfNN];
        while( neighbor != -1 && numberOfNN < maxNN){
            
            rhoTemp += m[neighbor]*cubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml);
            numberOfNN++;
            neighbor = NNsquare[maxNN * pCounter + numberOfNN];
            
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
        
        double gradKernel[DIM] = {0,0,0};
        
        double deltaVx;
        double deltaVy;
        double deltaVz;

        neighbor = NNsquare[ maxNN * pCounter + numberOfNN];
       
        while( neighbor != -1 && numberOfNN < maxNN){
            
            gradCubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml, gradKernel);
            deltaVx = vx[pCounter]- vx[neighbor];
            deltaVy = vy[pCounter]- vy[neighbor];
            deltaVz = vz[pCounter]- vz[neighbor];
            
            drhoTemp += m[neighbor]*(deltaVx*gradKernel[0] + deltaVy*gradKernel[1] + deltaVz*gradKernel[2])/rho[neighbor];
            ++numberOfNN;
            
            if(numberOfNN < maxNN){
                neighbor = NNsquare[ maxNN * pCounter + numberOfNN];
            }
            else{
                neighbor = -1;
            }
        }
        
        drho[pCounter] = rho[pCounter]*drhoTemp;     
      
               
    } 
    
};

// Calc pressure over EOS, here isothermal, p dependent on sound velocity c_s
void Particles::compPressure(){
    bool negativePressure = false;
    for( int pCounter = 0; pCounter < N; pCounter++){
        p[pCounter] = c_s*c_s*(rho[pCounter] - rhoRel);
       
        if(std::isnan(p[pCounter]) ) {
                Logger(WARN) <<" >  pressure is nan, particle: " << pCounter;
            }
        if(p[pCounter] < 0) {
                negativePressure = true;
            }
        
    }
    /* if( negativePressure){
        Logger(WARN) <<" >  pressure is negative somewhere ";
    } */
    
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
                        stress[indexSigma(pCounter, i, j)] = -p[pCounter] + S_ij[indexMatrix(i,j)];
                        /* if( i == 0){
                        std::cout << p[pCounter] << " + " << S_ij[i*DIM+j] << " = " << stress[pCounter*DIM*DIM + i*DIM+ j] << std::endl;
                        } */
                    }
                    else{
                        stress[indexSigma(pCounter, i, j)] =  S_ij[indexMatrix(i,j)];
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
        double gradKernel[DIM] = {0,0,0};

        double vParticle[DIM]  = {vx[pCounter], vy[pCounter], vz[pCounter]};

        neighbor = NNsquare[ maxNN * pCounter + numberOfNN ];
        
        // loop over nearest neighbors of each particle to calc dvi/dxj
        while(neighbor != -1 && numberOfNN < maxNN){
            gradCubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml, gradKernel);
            double vNeighbor[DIM] = {vx[neighbor], vy[neighbor], vz[neighbor]};

            for( int i = 0; i < DIM; i++){
                for( int j = 0; j < DIM; j++){
                partialV[indexSigma(pCounter, i, j)] += (m[neighbor]/rho[neighbor]) *( vNeighbor[i]-vParticle[i]) * gradKernel[j];                    
                }    
            }

            // switch to  next neighbor
            numberOfNN++;
            neighbor = NNsquare[ maxNN * pCounter + numberOfNN ];

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
                epsilon[indexMatrix(i,j)] = 0.5*(partialV[indexSigma(pCounter, i, j)]+partialV[indexSigma(pCounter, j, i)] );
                omega[indexMatrix(i,j)] = 0.5*(partialV[indexSigma(pCounter, i, j)]- partialV[indexSigma(pCounter, j, i)] );
            }
        }

        dS11[pCounter] = 2*mu* (epsilon[indexMatrix(0,0)]- (1/3)* epsilon[indexMatrix(0,0)]);
        dS12[pCounter] = 2*mu* epsilon[indexMatrix(0,1)];
        dS22[pCounter] = 2*mu* (epsilon[indexMatrix(1,1)]- (1/3)* epsilon[indexMatrix(1,1)]);
        dS13[pCounter] = 2*mu* epsilon[indexMatrix(0,2)];
        dS23[pCounter] = 2*mu* epsilon[indexMatrix(1,2)];

        for(int k = 0; k < DIM; k++){
            
            dS11[pCounter] += S_ij[indexMatrix(0,k)]*omega[indexMatrix(0,k)] + omega[indexMatrix(0,k)]* S_ij[indexMatrix(k,0)];
            dS12[pCounter] += S_ij[indexMatrix(0,k)]*omega[indexMatrix(1,k)] + omega[indexMatrix(0,k)]* S_ij[indexMatrix(k,1)];
            dS22[pCounter] += S_ij[indexMatrix(1,k)]*omega[indexMatrix(1,k)] + omega[indexMatrix(1,k)]* S_ij[indexMatrix(k,1)];
            dS13[pCounter] += S_ij[indexMatrix(0,k)]*omega[indexMatrix(2,k)] + omega[indexMatrix(0,k)]* S_ij[indexMatrix(k,2)];
            dS23[pCounter] += S_ij[indexMatrix(1,k)]*omega[indexMatrix(2,k)] + omega[indexMatrix(1,k)]* S_ij[indexMatrix(k,2)];
        }

    }

};

#endif

#if ARTIFICIAL_VISCOSITY 
    double Particles::compArtificialVisc( int pCounter, int neighbor){
        
        double deltaRx;
        double deltaRy;
        double deltaRz;
        double deltaVx;
        double deltaVy;
        double deltaVz;

        deltaRx = x[pCounter]- x[neighbor];
        deltaRy = y[pCounter]- y[neighbor];
        deltaRz = z[pCounter]- z[neighbor];

        deltaVx = vx[pCounter]- vx[neighbor];
        deltaVy = vy[pCounter]- vy[neighbor];
        deltaVz = vz[pCounter]- vz[neighbor];

        double r_ij = (deltaRx * deltaRx + deltaRy * deltaRy + deltaRz * deltaRz);
        double vr_ij = (deltaRx * deltaVx + deltaRy * deltaVy + deltaRz * deltaVz); // v_ij*r_ij
        double P_ij = 0;

        if(vr_ij >= 0){
            return P_ij;
        }
        else{  
            double mu_ij = sml* vr_ij/(r_ij+ epsilon*sml*sml);   

            if(std::isinf(mu_ij)){
                Logger(WARN) << "mu_ij is +- inf";
            } 
            P_ij = (-alpha*c_s*mu_ij+ beta* mu_ij*mu_ij) / (0.5*(rho[pCounter]+rho[neighbor]));
            return P_ij;
        }


    };
#endif

#if ARTIFICIAL_STRESS
    void Particles::compRij(int pCounter, int neighbor, double* R_ab ){
#if SOLIDS
        // R_ab = R_a + R_b
        double R_a[DIM*DIM] = {0,0,0,0,0,0,0,0,0};
        double R_b[DIM*DIM] = {0,0,0,0,0,0,0,0,0};

        for(int i = 0; i < DIM; i++){
            for(int j = 0; j < DIM; j++){
                if( stress[indexSigma(pCounter, i,j)] > 0){
                    R_a[indexMatrix(i,j)] = -epsilonAS * stress[indexSigma(pCounter, i,j)] /(rho[pCounter]*rho[pCounter]);
                }
                if( stress[indexSigma(neighbor, i,j)] > 0){
                    R_b[indexMatrix(i,j)] = -epsilonAS * stress[indexSigma(neighbor, i,j)] /(rho[neighbor]*rho[neighbor]);
                }
                
                R_ab[indexMatrix(i,j)] = R_a[indexMatrix(i,j)] + R_b[indexMatrix(i,j)];

            }
        }
#else
        double R_a = 0;
        double R_b = 0;

        if(p[pCounter]< 0){
            R_a = epsilonAS*abs(p[pCounter])/ (rho[pCounter]*rho[pCounter]);
        }
        if(p[neighbor]< 0){
            R_b = epsilonAS*abs(p[neighbor])/ (rho[neighbor]*rho[neighbor]);
        }
        R_ab[0] = R_a + R_b;

#endif

    };

    double Particles::compFn(int pCounter, int neighbor){
        double fn = 0;
        fn = (cubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml)/(cubicSpline(deltaP, sml)));
        fn = pow(fn,n);
        return fn;  

    };
#endif

// Calc acceleration through differential equation, iterate over nearest neighbors
// TO DO: self-acceleration? derivative for kernel zero for r_ij = 0. --> not necessary
void Particles::compAcceleration(){
    // iterate over all particles
    
     for( int pCounter =0; pCounter < N; pCounter++){

        bool artiViscIsnan = false;
        int numberOfNN = 0;
        int neighbor = 0;

        double axTemp = 0;
        double ayTemp = 0;
        double azTemp = 0;

        double gradKernel[DIM] = {0,0,0};

#if SOLIDS
        double dV[DIM] = {0,0,0};
#else
        double prefactor = 0;
#endif

        neighbor = NNsquare[ maxNN * pCounter + numberOfNN];


        // sum over nearest neighbors (while- loop)
        while( neighbor != -1 && numberOfNN < maxNN){
       

#if ARTIFICIAL_VISCOSITY
            double Pi_ij = compArtificialVisc(pCounter, neighbor);

            if(std::isnan(Pi_ij)){
                Logger(WARN) << "Artificial viscosity is nan";
                artiViscIsnan = true;
                
            }
#endif

#if ARTIFICIAL_STRESS
        double fn = compFn(pCounter, neighbor); // needed for Solids and gases/liquids

#endif

#if SOLIDS
            // has to be set to zero for every neighbor
            dV[0] = 0;
            dV[1] = 0;
            dV[2] = 0;
            
            gradCubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml, gradKernel);

#if ARTIFICIAL_STRESS
            double R_ij[DIM*DIM] = {0,0,0,0,0,0,0,0,0};
            compRij(pCounter, neighbor, R_ij);
#endif

            // calc dv_i for one neighbor with sum over stress tensor
            for( int i = 0; i < DIM; i++){
                
                for( int j = 0; j < DIM; j++){
                    dV[i] += (stress[indexSigma(pCounter, i, j)]/(rho[pCounter]*rho[pCounter]) + 
                    stress[indexSigma(neighbor, i, j)]/( rho[neighbor]*rho[neighbor]))*gradKernel[j];

#if ARTIFICIAL_STRESS
                    dV[i] += R_ij[indexMatrix(i,j)]*fn * gradKernel[j];
#endif

                    if(std::isnan(dV[i])){
                        Logger(WARN) << " dV for i,j nan" << i << ", "<< j;
                    }

#if ARTIFICIAL_VISCOSITY
                    // add artificial viscosity
                    if(i == j){
                        dV[i] += -Pi_ij*gradKernel[j];

                        // check for NaNs
                        if(std::isnan(Pi_ij)){
                            Logger(WARN) << "Pi_ij nan: " << Pi_ij;
                        }
                        if(std::isnan(gradKernel[j])){
                            Logger(WARN) << "Kernel nan: " << gradKernel[j];
                        }
                        if(std::isnan(gradKernel[j]*Pi_ij)){
                            Logger(WARN) << "Kernel* Pi_ij nan: " << gradKernel[j] << " P_ij" <<  Pi_ij;
                        }
                        if(std::isnan(dV[i])){
                        Logger(WARN) << " dV for i,j nan with artificial viscosity" << i << ", " << j;
                        }
                 }
#endif// ARTIFICIAL VISCOSITY
                   //  dV[i] += (stress[indexSigma(pCounter, i, j)]+stress[indexSigma(neighbor, i, j)])/(rho[pCounter]*rho[neighbor]) 
                    // *gradKernel[j];  
                    
                }
            }

            // add up components of the sum to the total sum for each nearest neighbor
            axTemp += m[neighbor]* dV[0];
            ayTemp += m[neighbor]* dV[1];
            azTemp += m[neighbor]* dV[2];


#else //SOLIDS == 0
            //1st Version of SPH equation
            //prefactor = - m[neighbor] *(p[neighbor]+p[pCounter])/(rho[pCounter]*rho[neighbor]);

            // 2nd Version of SPH equation
            prefactor = -m[neighbor]* (p[neighbor]/(rho[neighbor]*rho[neighbor]) + p[pCounter]/(rho[pCounter]*rho[pCounter]));

#if ARTIFICIAL_VISCOSITY
            prefactor += -m[neighbor]*Pi_ij;
#endif // ARTIFICIAL_VISCOSITY

#if ARTIFICIAL_STRESS
            double R_ij[1] = {0};
            compRij(pCounter, neighbor, R_ij);
            prefactor += -m[neighbor]* R_ij[0]*fn;
#endif
            
            gradCubicSpline(x[pCounter], y[pCounter], z[pCounter], x[neighbor], y[neighbor], z[neighbor], sml, gradKernel);
            
            axTemp += prefactor*gradKernel[0];
            ayTemp += prefactor*gradKernel[1];
            azTemp += prefactor*gradKernel[2];
#endif // SOLIDS

            numberOfNN++;
            neighbor = NNsquare[ maxNN * pCounter + numberOfNN]; // get next neighbors

        } // end of while-loop

        if(artiViscIsnan){
            Logger(WARN) << "artificial viscosity nan somewhere... ";
        }
        if(std::isnan(axTemp) || std::isnan(ayTemp) || std::isnan(azTemp) ) {
                Logger(WARN) <<" >  acceleration is nan, particle: " << pCounter;
                //exit(1);
                //printf("ax: %f, ay: %f, az: %f \n", axTemp, ayTemp ,azTemp);
                //printf(" artificial vis: %f :" , P_ij)
            }

        // write new accelerations 
        ax[pCounter] = axTemp;
        ay[pCounter] = ayTemp;
        az[pCounter] = azTemp;       
               
    } // end of for-loop
 
};

#if COURANT_CONDITION
    double Particles::compTimestep(double timestep){
        double maxStep = 0;
#if ARTIFICIAL_VISCOSITY

        double mu_max = std::numeric_limits<double>::min(); // smallest possible value for mu_max

        for(int pCounter =0; pCounter < N; pCounter++){
            int numberOfNN = 0;
            int neighbor = NNsquare[ pCounter * maxNN + numberOfNN];

            while(neighbor != -1 && neighbor < maxNN){

                double deltaRx;
                double deltaRy;
                double deltaRz;
                double deltaVx;
                double deltaVy;
                double deltaVz;

                deltaRx = x[pCounter]- x[neighbor];
                deltaRy = y[pCounter]- y[neighbor];
                deltaRz = z[pCounter]- z[neighbor];

                deltaVx = vx[pCounter]- vx[neighbor];
                deltaVy = vy[pCounter]- vy[neighbor];
                deltaVz = vz[pCounter]- vz[neighbor];

                double r_ij = (deltaRx * deltaRx + deltaRy * deltaRy + deltaRz * deltaRz);
                double vr_ij = (deltaRx * deltaVx + deltaRy * deltaVy + deltaRz * deltaVz);

                double mu = sml* vr_ij/(r_ij+ epsilon*sml*sml);   
                mu_max = mu_max < mu ? mu : mu_max;

                // next neighbor
                numberOfNN++;
                neighbor = NNsquare[pCounter *maxNN + numberOfNN];

            }
        }
        maxStep = CFL*sml/(c_s+ 1.2*(c_s*alpha + mu_max*beta));

#else // no Artificial viscosity
        maxStep = CFL*sml/c_s;

#endif //ARTIFICIAL_VISCOSITY
        if(maxStep < timestep){
            return maxStep;
        }
        else{
            return timestep;
        }
        
    }

#endif // COURANT_CONDITION

void Particles::write2file(std::string filename){
    // open output file
    HighFive::File h5File { filename, HighFive::File::ReadWrite |
                                      HighFive::File::Create |
                                      HighFive::File::Truncate };
    
    // dimensions for datasets containing vectors
    std::vector<size_t> dataSpaceDims(2);
    dataSpaceDims[0] = std::size_t(N); // number of particles
    dataSpaceDims[1] = DIM;

#if SOLIDS
    // dimension of tensor
    std::vector<size_t> dataSpaceDims2(2);
    dataSpaceDims2[0] = std::size_t(N); // number of particles
    dataSpaceDims2[1] = DIM*DIM;

#endif

#if NN
    std::vector<size_t> dataSpaceDims3(2); // for nearest neighbors
    dataSpaceDims3[0] = std::size_t(N);
    dataSpaceDims3[1] = maxNN;
#endif


    // create datasets
    HighFive::DataSet rhoDataSet = h5File.createDataSet<double>("/rho", HighFive::DataSpace(N));
    HighFive::DataSet pDataSet = h5File.createDataSet<double>("/p", HighFive::DataSpace(N));
    HighFive::DataSet mDataSet = h5File.createDataSet<double>("/m", HighFive::DataSpace(N));
    HighFive::DataSet iCounterDataSet = h5File.createDataSet<int>("/iCounter", HighFive::DataSpace(N)); // count of interaction partners
    
#if CALC_DENSITY
    HighFive::DataSet drhoDataSet = h5File.createDataSet<double>("/drho", HighFive::DataSpace(N));
#endif

#if SOLIDS
    HighFive::DataSet SDataSet = h5File.createDataSet<double>("/S", HighFive::DataSpace(dataSpaceDims2)); // deviatoric stress   
    HighFive::DataSet stressDataSet = h5File.createDataSet<double>("/sigma", HighFive::DataSpace(dataSpaceDims2)); // stress tensor
    HighFive::DataSet dvDataSet = h5File.createDataSet<double>("/dvdx", HighFive::DataSpace(dataSpaceDims2));  // partial derivatives dv/dx
#endif
    
#if NN
    HighFive::DataSet NNDataSet = h5File.createDataSet<int>("/NN", HighFive::DataSpace(dataSpaceDims3));
#endif

    HighFive::DataSet posDataSet = h5File.createDataSet<double>("/r", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet velDataSet = h5File.createDataSet<double>("/v", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet accDataSet = h5File.createDataSet<double>("/a", HighFive::DataSpace(dataSpaceDims));

      // containers for particle data
    std::vector<double> rhoVec(rho,rho+N);
    std::vector<double> pVec(p,p+N);
    std::vector<double> mVec(m, m + N);

    std::vector<int> iCounterVec(iCounter, iCounter+N);

#if CALC_DENSITY
    std::vector<double> drhoVec(rho,rho+N);
#endif
    

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

#if NN
    std::vector<std::vector<int>> NNVec(N);
    std::vector<int> NNbuf(maxNN);
    for( int pCounter = 0; pCounter < N; pCounter++){
        for(int nCounter = 0; nCounter < maxNN; nCounter++){
            NNbuf[nCounter] = NNsquare[ maxNN * pCounter + nCounter];
        }
        NNVec[pCounter] = NNbuf;
    }
#endif

#if SOLIDS
    std::vector<std::vector<double>> stressVec(N);
    std::vector<std::vector<double>> dvVec(N);
    std::vector<std::vector<double>> sVec(N);

    std::vector<double> stressBuf(DIM*DIM);
    std::vector<double> dvBuf(DIM*DIM);
    std::vector<double> sBuf(DIM*DIM);

    for( int pCounter = 0; pCounter < N; pCounter++){
        for(int i = 0; i < DIM; i++){
            for( int j = 0; j < DIM; j++){
                stressBuf[i*DIM+j] = stress[pCounter*DIM*DIM + i*DIM + j];
                dvBuf[i*DIM+j] = partialV[pCounter*DIM*DIM + i*DIM + j];
            }
        }
        sBuf[0] = S11[pCounter];
        sBuf[1] = S12[pCounter];
        sBuf[2] = S13[pCounter];
        sBuf[3] = S12[pCounter];
        sBuf[4] = S22[pCounter];
        sBuf[5] = S23[pCounter];
        sBuf[6] = S13[pCounter];
        sBuf[7] = S23[pCounter];
        sBuf[8] = -(S11[pCounter]+ S22[pCounter]);

        stressVec[pCounter] = stressBuf;
        dvVec[pCounter] = dvBuf;
        sVec[pCounter] = sBuf;
       
    }
#endif
    // write data
    rhoDataSet.write(rhoVec);
    pDataSet.write(pVec);
    mDataSet.write(mVec);
    iCounterDataSet.write(iCounterVec);

#if CALC_DENSITY
    drhoDataSet.write(drho);
#endif
    
    posDataSet.write(posVec);
    velDataSet.write(velVec);
    accDataSet.write(accVec);  
#if NN
    NNDataSet.write(NNVec);
#endif 

#if SOLIDS
    SDataSet.write(sVec);
    dvDataSet.write(dvVec);
    stressDataSet.write(stressVec);
#endif

};

