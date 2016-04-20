

#ifndef INTERPOLATE_H_INCLUDED
#define INTERPOLATE_H_INCLUDED

#include <stdlib.h>
#include "frame.h"
#include "fd_der.h"
#define FD_MODE 0		// 0 MEANS JUST SIMPLE CENTRAL 2ND ORDER FINITE DIFFERENCE FOR 2ND DERIVATIVE
				// 1 MEANS FANCY THREE TYPES OF DERIVATIVES AND PICK 

extern int N;
extern int SIZE;
extern double H;



double cal_interpgrad(double xj, double yj, double **gradx, double **grady, double *interpolatex, double *interpolatey);

double cal_curv(double xi, double xj, double **dist);

double cal_dist(double x, double y, double **dist, double **gradx, double **grady);


#endif
