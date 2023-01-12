#include <cmath>
#include "Kernel.h"

#define _USE_MATH_DEFINE

// see NumHydro script

// calc cubic Spline Kernel
double cubicSpline(double xi, double yi, double zi, double xj, double yj, double zj, double sml){
    // one h 
    double rij= sqrt((xi-xj)*(xi-xj) + (yi-yj)*(yi-yj)+ (zi-zj)*(zi-zj));
    double rh = rij/sml;

    double prefactor = 8/(M_PI*pow(sml,DIM)); //depends on DIM
    
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
double* gradCubicSpline(double xi, double yi, double zi, double xj, double yj, double zj, double sml){
    // one h
    double xij = xi - xj;
    double yij = yi - yj;
    double zij = zi - zj;
    double rij = sqrt(xij*xij + yij*yij + zij*zij);
    double rh = rij/sml;
    
    double gradW;
    double prefactor = 6*8/(M_PI*pow(sml,DIM+1));

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

    static double gradWvec[DIM]; //static necessary

    gradWvec[0] = gradW* xij/rij;
    gradWvec[1] = gradW* yij/rij;
    gradWvec[2] = gradW* zij/rij;
    
    return gradWvec;
}