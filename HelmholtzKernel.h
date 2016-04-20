
#ifndef HELMHOLTZ_KERNEL_H_DEFINED
#define HELMHOLTZ_KERNEL_H_DEFINED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "constants.h"
#include "frame.h"
#include "tools.h"
#include "functions.h"
#include "PDE.h"
#include "fd_der.h"
#include "interpolate.h"

extern int N;
extern double H, HSQUARED;
extern int SIZE;
extern double WAVE;

int new_HelmholtzKernel(int size, double delta, double epsilon, \
	double *zxlist, double *zylist, double *zxstarlist, \
	double *zystarlist, double *dist, double *gradx, double *grady, double *curv, double *Jacobian, \
	double *PolyW, double *SineW, double *PolyResult, double *SineResult, \
	double **PolyKernelH, double **PolyKernelHimg, double **SineKernelH, double **SineKernelHimg);

int new_HelmholtzKernelCombo(int size, double delta, double epsilon, double **level, double **leveldx, double **leveldy, \
	double **levelcurv, double *zxlist, double *zylist, double *zxstarlist, \
	double *zystarlist, double *dist, double *gradx, double *grady, double *curv, double *Jacobian, double *regfactor, \
	double *PolyW, double *SineW, double *PolyResult, double *SineResult, \
	double **PolyKernelH, double **PolyKernelHimg, double **SineKernelH, double **SineKernelHimg, double wavenum, \
	double eta);


#endif
