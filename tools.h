
#ifndef TOOLS_H_DEFINED
#define TOOLS_H_DEFINED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "constants.h"

#define BESSELJ_APROX 4.0	//	determine the number of terms for approximating J0	//
#define BESSELY_APROX 4.0	//	determine the number of terms for approximating Y0	//
				//	BOTH OF THE ABOVE DEPEND ON THE WAVE NUMBER	//

#define ROOTBESSEL {{0.0, 2.404825558, 5.520078110, 8.653727913, 11.791534439, 14.930917708}, \
		    {0.0, 3.831705970, 7.015586670, 10.17346814, 13.323691936, 16.470630051}, \
		    {0.0, 5.135622302, 8.417244140, 11.61984117, 14.795951782, 17.959819495}, \
		    {0.0, 6.380161896, 9.761023130, 13.01520072, 16.223466160, 19.409415226}, \
		    {0.0, 7.588342435, 11.06470949, 14.37253667, 17.615966050, 20.826932957}, \
		    {0.0, 8.771483816, 12.33860420, 15.70017408, 18.980133875, 22.217799897} }


//	for BESSEL FUNCTIONS 0 < x < 3	//
#define J0COEFFICIENT {1.0, -2.249999879, 1.26562306, -0.316394552, 0.044460948, -0.003954479, 0.000212950}
#define Y0COEFFICIENT {0.367466907, 0.605593797, -0.743505078, 0.253005481, -0.042619616, 0.004285691, -0.000250716}
#define J1COEFFICIENT {0.5, -0.5625, 0.210937377, -0.03955004, 0.004447331, -0.000330547, 0.000015525}
#define Y1COEFFICIENT {0.073735531, 0.722769344, -0.438896337, 0.104320251, -0.013637596, 0.00112597, -0.000056455}

//	for BESSEL FUNCTIONS x > 3	//
#define F0COEFFICIENT {0.79788454, -0.00553897, 0.00099336, -0.00044346, 0.00020445, -0.0000496}
#define THETA0COEFFICIENT {-0.04166592, 0.002394, -0.00073984, 0.000311, -0.00007605, 0.0}
#define F1COEFFICIENT {0.79788459, 0.01662008, -0.00187, 0.0006852, -0.0002944, 0.00006952}
#define THETA1COEFFICIENT {0.12499895, -0.0060524, 0.00135825, -0.00049616, 0.00011531, 0.0}


double inner(int size, double vectorA[], double vectorB[]);
double omp_inner(int size, double vectorA[], double vectorB[], int chunk);
double norm(int size, double vectorA[]);
double min(double a, double b);
double max(double a, double b);
double minmod(double x, double y);
double min_abs(double a, double b, double c);
double Godunov(double signage, double a, double b, double c, double d);

double BesselJ0(double x);
double BesselY0(double x);
double BesselJ1(double x);
double BesselY1(double x);
double BesselJ(int order, double x);
double BesselY(int order, double x);
double besselj(int order, double x);
double bessely(int order, double x);
double factorial(double n);




#endif
