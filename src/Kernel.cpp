//
// created by Anne Vera Jeschke December 2022
//
#include <cmath>
#include "Kernel.h"

#define _USE_MATH_DEFINE

// see NumHydro script

// calc cubic Spline Kernel
double cubicSpline(double xi, double yi, double zi, double xj, double yj, double zj, double sml){
    // one h 
    double rij= sqrt((xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj)  );
    double rh = rij/sml;

    double prefactor = 0.0;

#if KERNEL_DIM == 3
    prefactor = 8/(M_PI*pow(sml,KERNEL_DIM)); //depends on KERNEL_DIM
#else //KERNEL_DIM == 2
    prefactor = 40/(7*M_PI*sml*sml);
#endif

    if(rh > 1){
        return 0.0;
    }
    else if (rh >= 0.5 && rh <= 1)
    {
        return 2*pow((1-rh), 3)*prefactor;
    }
    else{
        return (6*pow(rh,3)-6*rh*rh +1)*prefactor;
    }
    
    
}

double cubicSpline(double r, double sml){
    double rij= r;
    double rh = rij/sml;

    double prefactor = 0.0;

#if KERNEL_DIM == 3
    prefactor = 8/(M_PI*pow(sml,KERNEL_DIM)); //depends on KERNEL_DIM
#else //KERNEL_DIM == 2
    prefactor = 40/(7*M_PI*sml*sml);
#endif

    if(rh > 1){
        return 0.0;
    }
    else if (rh >= 0.5 && rh <= 1)
    {
        return 2*pow((1-rh), 3)*prefactor;
    }
    else{
        return (6*pow(rh,3)-6*rh*rh +1)*prefactor;
    }

}

// calc gradient of Cubic Spline Kernel
void gradCubicSpline(double xi, double yi, double zi, double xj, double yj, double zj, double sml, double *grad){
    // one h
    double xij = xi - xj;
    double yij = yi - yj;
    double zij = zi - zj;

    double rij = sqrt(xij*xij + yij*yij + zij*zij );

    double rh = rij/sml;
    
    double gradW;
#if KERNEL_DIM == 3
    double prefactor = 6*8/(M_PI*pow(sml,KERNEL_DIM+1));
#else // KERNEL_DIM == 2
    double prefactor = 6*40/(7*M_PI*pow(sml, KERNEL_DIM+1));
#endif

    if(rh > 1){
        gradW = 0.0;
    }
    else if (rh >= 0.5 && rh <= 1)
    {
        gradW = -(1-rh)*(1-rh);
    }
    else{
        gradW = 3*rh*rh-2*rh;
    }
    gradW = gradW*prefactor;
  
    grad[0] = gradW* xij/rij;
    grad[1] = gradW* yij/rij;
    grad[2] = gradW* zij/rij;

    
}