
#include "FarFieldEval.h"

double FarField(double x1, double x2, double delta0, double epsilon0, double *u, double *uimg)
{
	int i, j;
	double pointxnorm, pointynorm, directionx1, directionx2, directiony1, directiony2;
//	double
	double constterm;
	double tempx, tempy, tstarx, tstary;

	pointxnorm = sqrt(x1*x1 + x2*x2);
	directionx1 = x1/pointxnorm;
	directionx2 = x2/pointxnorm;

	consttermR = 0.25/sqrt(PI*WAVE);	//	{ e^(ipi/4) }/sqrt(k pi)	//


	for (i = 0; i < SIZE; i++)
	{
		tempy = INT_L + (double)(i) * H;
	for (j = 0; j < SIZE; j++)
	{
		tempx = INT_L + (double)(j) * H;
		
		thisdist = DIST[i][j];
		if (TRAD_ONESIDED == 0)
		{
			if (fabs(thisdist) > epsilon0)
				continue;
			delta = cal_delta(epsilon0, thisdist);
		}
		else
		{
			if ( (thisdist > epsilon0) || (thisdist < delta0))
				continue;
			delta = 2. * cal_delta(epsilon0, thisdist);
		}
		tstarx = tempx - DIST_DX[i][j] * thisdist;
		tstary = tempy - DIST_DY[i][j] * thisdist;

		pointynorm = sqrt(tstarx*tstarx + tstary*tstary);
		directiony1 = tstarx/pointynorm;
		directiony2 = tstary/pointynorm;

		cal_interpgrad(tstarx, tstary, DIST_DX, DIST_DY, &stargradx, &stargrady);
		stargradnorm = sqrt(stargradx*stargradx + stargrady*stargrady);
		stargradx = stargradx/stargradnorm;
		stargrady = stargrady/stargradnorm;

		thisJ = 1. + thisdist * CURV[i][j];
		result = H*H*delta*thisJ;

		inner_pd1 = -1. * (directionx1 * stargradx + directionx2 * stargrady);
		inner_pd2 = (directionx1 *directiony1 + directionx2 * directiony2);

		KernelR = ( WAVE * inner_pd1 )
	
		

	}
	}

}
