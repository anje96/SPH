class Particles{
    public: 
        Particles(int NumParticles);
        ~Particles();

        int N;
        double *m, *x, *y, *z, *vx, *vy, *vz, *rho, *p, *ax, *ay, *az;
        
};

