
#ifndef FD_DER_H_INCLUDED
#define FD_DER_H_INCLUDED

#include <stdlib.h>
#include "tools.h"

extern int N, SIZE;
extern double H;



double XDER_C(double **level, int i, int j, int mode);
double YDER_C(double **level, int i, int j, int mode);

double XDER2(double **level, int i, int j, int mode);
double YDER2(double **level, int i, int j, int mode);

void FD_DX(double **phi, int i, int j, double *Dxp, double *Dxn, int mode);
void FD_DY(double **phi, int i, int j, double *Dyp, double *Dyn, int mode);







#endif
