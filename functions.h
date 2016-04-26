
#ifndef FUNCTIONS_H_DEFINED
#define FUNCTIoNS_H_DEFINED


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//	#include <gsl_sf_bessel.h>

#include "PDE.h"
#include "frame.h"
#include "constants.h"
#include "tools.h"

			//	BOTH OF THE ABOVE DEPEND ON THE WAVE NUMBER	//
#define PSI_ETA 0.012	//	determine how smooth the point source is, the smaller the more concentrated 	//
#define DIST_MODE 0 	//	0 for circle(s), 1 for level set function	//
#define NONHOM_MODE 0	//	1 for nonhomogeneous mode, 0 for homogeneous mode.	//
#define MOMENT_ORDER 1
#define AVG_KERNEL "cosine"

#define CONSTA1 -759.2781934172483
#define CONSTB1 446.2604260472818

#define CONSTA2 14317.969703708994
#define CONSTB2 -16509.044867497203
#define CONSTC2 4480.224717878304

extern int N;
extern double H;
extern int SIZE;
extern double WAVE;

//	double cal_J(int mode, double x, double y);
double cal_partial(int mode, double x1, double y1, double x2, double y2);
double cal_regpartial(double x1, double y1, double x2, double y2, double tau);
double cal_delta(double epsilon, double x);
//	double cal_interpgrad(int mode, double xj, double yj);
double cal_nonhom(double zstarx, double zstary, double A[]);

double gradientdx(double x, double y);
double gradientdy(double x, double y);
double laplace(double x, double y);
double phi(double xi, double yi, double xj, double yj);
double phiimg(double xi, double yi, double xj, double yj);
double psi(double xj, double yj);
double psiimg(double xj, double yj);


//	double starx(double zx, double zy);
//	double stary(double zx, double zy);

double distance(double x, double y);

double PolyWeight(double t, double tau);
double SineWeight(double t, double tau);

#endif
