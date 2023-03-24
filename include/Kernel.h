#ifndef SPH_KERNEL_H
#define SPH_KERNEL_H

#include "parameter.h"


double cubicSpline(double xi, double yi, double zi, double xj, double yj, double zj, double sml);

void gradCubicSpline(double xi, double yi, double zi, double xj, double yj, double zj, double sml, double *grad);

#endif // SPH_KERNEL_H