
#include "HelmholtzKernel.h"

int new_HelmholtzKernel(int size, double delta, double epsilon, double *zxlist, double *zylist, double *zxstarlist, \
	double *zystarlist, double *dist, double *gradx, double *grady, double *curv, double *Jacobian, \
	double *PolyW, double *SineW, double *PolyResult, double *SineResult, \
	double **PolyKernelH, double **PolyKernelHimg, double **SineKernelH, double **SineKernelHimg)
{
	int i, j;
	int counter;
	double zx, zy, zstarx, zstary, thisdist, thisgradx, thisgrady, stargradx, stargrady;
	double xminusy1, xminusy2;
	double starxi, starxj, staryi, staryj, partialx, partialy, partialximg, partialyimg, partial, partialimg;

	double regfactor, regfactorimg, threshold, tau;

	threshold = 0.01*H*H;

	tau = H/TAU_FACTOR;
	threshold = tau*tau;

	counter = 0;
	for (i = 0; i < SIZE; i++)
	{
		zy = INT_L + (double)(i) * H;
		for (j = 0; j < SIZE; j++)
		{	
			zx = INT_L + (double)(j) * H;

			thisdist = radius - sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) );
			
			if (COMBO_ONESIDED == 0)
			{
				if ( ( fabs(thisdist) > epsilon) || (fabs(thisdist) < delta) )
					continue;
			}
			else
			{
				if ( ( thisdist > epsilon) || (thisdist < delta) )
					continue;
			}

			thisgradx = (Cx - zx) / sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) );
			thisgrady = (Cy - zy) / sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) );
			
			zstarx = zx - thisdist * thisgradx;
			zstary = zy - thisdist * thisgrady;

			stargradx = (Cx - zstarx) / sqrt( (zstarx-Cx)*(zstarx-Cx) + (zstary-Cy)*(zstary-Cy) );
			stargrady = (Cy - zstary) / sqrt( (zstarx-Cx)*(zstarx-Cx) + (zstary-Cy)*(zstary-Cy) );

			zxlist[counter] = zx;
			zylist[counter] = zy;
			zxstarlist[counter] = zstarx;
			zystarlist[counter] = zstary;
			dist[counter] = thisdist;
			gradx[counter] = stargradx;
			grady[counter] = stargrady;
			curv[counter] = 1./(radius - thisdist);
			
			Jacobian[counter] = 1. + thisdist * curv[counter];
//			Jacobian[counter] = 1.;

			if (POLY_TEST == 1)
			{
				if (COMBO_ONESIDED == 0)
					PolyW[counter] = PolyWeight(fabs(thisdist/epsilon), delta/epsilon)/(2.0*epsilon);
				else
					PolyW[counter] = PolyWeight(fabs(thisdist/epsilon), delta/epsilon)/(epsilon);
				PolyResult[counter] = H*H*Jacobian[counter] * PolyW[counter];
			}
			if (SINE_TEST == 1)
			{
				if (COMBO_ONESIDED == 0)
					SineW[counter] = SineWeight(fabs(thisdist/epsilon), delta/epsilon)/(2.0*epsilon);
				else
					SineW[counter] = SineWeight(fabs(thisdist/epsilon), delta/epsilon)/(epsilon);
				SineResult[counter] = H*H*Jacobian[counter] * SineW[counter];
			}

			counter++;		
		}
	}

	if (counter != size)
	{
		printf("Size different for new Kernel.\n");
		exit(0);
	}

	for (i = 0; i < size; i++)
	{
		starxi = zxstarlist[i];
		staryi = zystarlist[i];
		for ( j = 0; j < size; j++)
		{
			starxj = zxstarlist[j];
			staryj = zystarlist[j];

			partialx = cal_partial(1, starxi, staryi, starxj, staryj);// x coordinate real	//
			partialy = cal_partial(2, starxi, staryi, starxj, staryj);// y coordinate real	//
			partialximg = cal_partial(3, starxi, staryi, starxj, staryj);// x coordinate imaginary //
			partialyimg = cal_partial(4, starxi, staryi, starxj, staryj);// y coordinate imaginary //

			partial =  -1. * (partialx * gradx[i] + partialy * grady[i]);
			partialimg = -1. * (partialximg * gradx[i] + partialyimg * grady[i]);

			if ( ( (starxi-starxj)*(starxi-starxj) + (staryi-staryj)*(staryi-staryj) )< threshold )
			{
				partial = -1.*curv[j]/(4.*PI);
				partialimg = 0.;
			}
//			else
//			{
//				partial = partial;
//				partialimg = partialimg;
//			}
			if (POLY_TEST == 1)
			{
				PolyKernelH[i][j] = partial * PolyResult[j];
				PolyKernelHimg[i][j] = partialimg * PolyResult[j];
				if (i == j)
					PolyKernelH[i][j] -= 0.5;
			}
			if (SINE_TEST == 1)
			{
				SineKernelH[i][j] = partial * SineResult[j];
				SineKernelHimg[i][j] = partialimg * SineResult[j];
				if (i == j)
					SineKernelH[i][j] -= 0.5;
			}
		}
	}
	return 0;
}



int new_HelmholtzKernelCombo(int size, double delta, double epsilon, double **level, double **leveldx, double **leveldy, \
	double **levelcurv, double *zxlist, double *zylist, double *zxstarlist, \
	double *zystarlist, double *dist, double *gradx, double *grady, double *curv, double *Jacobian, double *regfactor, \
	double *PolyW, double *SineW, double *PolyResult, double *SineResult, \
	double **PolyKernelH, double **PolyKernelHimg, double **SineKernelH, double **SineKernelHimg, double wavenum, \
	double eta)
{
	int i, j;
	int counter;
	double zx, zy, zstarx, zstary, thisdist, thisgradx, thisgrady, stargradx, stargrady, stargradnorm;
	static const double recipfourpi = 1./(4.*PI);

//	double eta;

	double threshold, tau, regterm1, regterm2, regterm3, regterm4;

	tau = H/TAU_FACTOR;
	threshold = tau * tau;


//	if (INTERIOR_MODE == 1)
//		regfactor = 1./(4.*PI*radius);
//	else
//		regfactor = -1./(4.*PI*radius);
//	regfactorimg = 0.0;

//	threshold = 0.01*h*h;

////////////////////////////////////////////////////////////
//	eta = 0.5;
////////////////////////////////////////////////////////////

	counter = 0;
	for ( i = 0; i < SIZE; i++)
	{
		zy = INT_L + (double)(i) * H;
		for (j = 0; j < SIZE; j++)
		{	
			zx = INT_L + (double)(j) * H;

//			thisdist = radius - sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) );

			thisdist = level[i][j];	

			if (COMBO_ONESIDED == 0)
			{		
				if ( ( fabs(thisdist) > epsilon) || (fabs(thisdist) < delta) )
					continue;
			}
			else
			{
				if ( ( thisdist > epsilon) || (thisdist < delta) )
				continue;
			}

//			thisgradx = (Cx - zx) / sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) );
//			thisgrady = (Cy - zy) / sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) );

			thisgradx = leveldx[i][j];
			thisgrady = leveldy[i][j];
			
			zstarx = zx - thisdist * thisgradx;
			zstary = zy - thisdist * thisgrady;

//			stargradx = (Cx - zstarx) / sqrt( (zstarx-Cx)*(zstarx-Cx) + (zstary-Cy)*(zstary-Cy) );
//			stargrady = (Cy - zstary) / sqrt( (zstarx-Cx)*(zstarx-Cx) + (zstary-Cy)*(zstary-Cy) );

			cal_interpgrad(zstarx, zstary, leveldx, leveldy, &stargradx, &stargrady);
//			printf("stargrad (%lf, %lf)\n", stargradx, stargrady);
			stargradnorm = sqrt(stargradx*stargradx + stargrady*stargrady);
			stargradx = stargradx/stargradnorm;
			stargrady = stargrady/stargradnorm;

//			if ( (stargradx != stargradx) || (stargrady != stargrady) )
//			{
//				printf("(%d, %d), point (%lf, %lf), dist %lf, grad (%lf, %lf)\n", \
//				i, j, zx, zy, thisdist, thisgradx, thisgrady);
//				printf("Projected to (%lf, %lf), grad (%lf, %lf), norm %lf\n", \
//				zstarx, zstary, stargradx, stargrady, stargradnorm);
//				exit(0);
//			}
			
			zxlist[counter] = zx;
			zylist[counter] = zy;
			zxstarlist[counter] = zstarx;
			zystarlist[counter] = zstary;
			dist[counter] = thisdist;
			gradx[counter] = stargradx;
			grady[counter] = stargrady;
//			curv[counter] = 1./(radius - thisdist);
			curv[counter] = cal_curv(zstarx, zstary, level);
			regfactor[counter] = curv[counter]*recipfourpi;
			
			Jacobian[counter] = 1. + thisdist * levelcurv[i][j];
//			Jacobian[counter] = 1.;

			if (POLY_TEST == 1)
			{
				if (COMBO_ONESIDED == 0)
					PolyW[counter] = PolyWeight(fabs(thisdist/epsilon), delta/epsilon)/(2.0*epsilon);
				else
					PolyW[counter] = PolyWeight(fabs(thisdist/epsilon), delta/epsilon)/(epsilon);
				PolyResult[counter] = HSQUARED*Jacobian[counter] * PolyW[counter];
			}
			if (SINE_TEST == 1)
			{
				if (COMBO_ONESIDED == 0)
					SineW[counter] = SineWeight(fabs(thisdist/epsilon), delta/epsilon)/(2.0*epsilon);
				else
					SineW[counter] = SineWeight(fabs(thisdist/epsilon), delta/epsilon)/(epsilon);
				SineResult[counter] = HSQUARED*Jacobian[counter] * SineW[counter];
			}
//			if (SineResult[counter] != SineResult[counter])
//			{
///				printf("(i, j) = (%d, %d), point (%lf, %lf), star (%lf ,%lf) dist %lf\n", \
///				i, j, zx, zy, zstarx, zstary, thisdist);
//				printf("grad (%lf, %lf), curv %lf, Weight %lf, J %lf\n", \
//				gradx[counter], grady[counter], curv[counter], SineW[counter], Jacobian[counter]);
//				exit(0);
//			}

			counter++;		
		}
	}

	if (counter != size)
	{
		printf("Size different for new Kernel.\n");
		exit(0);
	}

	#pragma omp parallel for private(i,j) schedule(static, 16)
	for (i = 0; i < size; i++)
	{
		double xminusy1, xminusy2, xminusynorm, recipxminusynorm;
		double starxi, starxj, staryi, staryj, partialx, partialy, partialximg, partialyimg, partial, partialimg;
		double Polysinglereal, Polysingleimg, Polydoublereal, Polydoubleimg, Polytempreal, Polytempimg;
		double Sinesinglereal, Sinesingleimg, Sinedoublereal, Sinedoubleimg, Sinetempreal, Sinetempimg;

		double x1, x2, parameter, BesselJ1value, BesselY1value, commonterm1, commonterm2;
		double inner_pd1, inner_pd2, inner_pd3;
		double term1, term2, term3, term1R, term2R, term3R, term1I, term2I, term3I;

		double xmy1, xmy2, xmynorm, recipxmynorm;
		double single_BesselJ1value, single_BesselY1value, single_commonterm, single_parameter;
		double xmyinner_pd;

		starxi = zxstarlist[i];
		staryi = zystarlist[i];
		for ( j = 0; j < size; j++)
		{
			
			starxj = zxstarlist[j];
			staryj = zystarlist[j];

			//	SINGLE LAYER PART, dG(x*,y*)/dnx W(y) J(y) * h^2	//

			xmy1 = starxi-starxj;
			xmy2 = staryi-staryj;
			xmynorm = sqrt(xmy1*xmy1 + xmy2*xmy2);

			if (xmynorm > tau)
			{
				recipxmynorm = 1./xmynorm;
				parameter = WAVE * xmynorm;
				single_commonterm = 0.25 * WAVE * recipxmynorm;
				single_BesselJ1value = BesselJ1(parameter);
				single_BesselY1value = BesselY1(parameter);

//				partialx = single_commonterm * single_BesselY1value * xmy1;
//				partialy = single_commonterm * single_BesselY1value * xmy2;
//				partialximg = -1. * single_commonterm * single_BesselJ1value * xmy1;
//				partialyimg = -1. * single_commonterm * single_BesselJ1value * xmy2;
//				partial =  -1. * (partialx * gradx[i] + partialy * grady[i]);
//				partialimg = -1. * (partialximg * gradx[i] + partialyimg * grady[i]);

				xmyinner_pd = (xmy1 * gradx[i] + xmy2 * grady[i]);
				partial = -1. * single_commonterm * single_BesselY1value * xmyinner_pd;
				partialimg = single_commonterm * single_BesselJ1value * xmyinner_pd;
				
			}
//			if (0)
//			{
//			partialx = cal_partial(1, starxi, staryi, starxj, staryj);// x coordinate real	//
//			partialy = cal_partial(2, starxi, staryi, starxj, staryj);// y coordinate real	//
//			partialximg = cal_partial(3, starxi, staryi, starxj, staryj);// x coordinate imaginary //
//			partialyimg = cal_partial(4, starxi, staryi, starxj, staryj);// y coordinate imaginary //
//			}
			//	ik/4 * H1(k|x-y|)/|x-y| [(x-y) dot nx] = ik/4 * (J1 + iY1)/|x-y| [inner_pd]	//
			else
			{
				partial = regfactor[i];
				partialimg = 0.0;
			}			

//			if ( ( (starxi-starxj)*(starxi-starxj) + (staryi-staryj)*(staryi-staryj) )< threshold )
//			{
//				partial = -1.*curv[j]/(4.*PI);
//				partialimg = 0.;
//			}
//			else
//			{
//				
//			}
//			else
//			{
//				partial = partial;
//				partialimg = partialimg;
//			}
			if (POLY_TEST == 1)
			{
				Polytempreal = partial * PolyResult[j];
				Polytempimg = partialimg * PolyResult[j];
				if (i == j)
				{
					Polytempreal -= 0.5;
//					Polytempreal -= 1.0;
				}
				//	-i * eta * dG(x,y*)/dnx		//
				Polysinglereal = eta * Polytempimg;
				Polysingleimg = -1. * eta * Polytempreal;
			}
			if (SINE_TEST == 1)
			{
				Sinetempreal = partial * SineResult[j];
				Sinetempimg = partialimg * SineResult[j];
				if (i == j)
				{
					Sinetempreal -= 0.5;
//					Sinetempreal -= 1.0;
				}
				Sinesinglereal = eta * Sinetempimg;
				Sinesingleimg = -1. * eta * Sinetempreal;
			}

			//	DOUBLE LAYER PART, d^2G(x(d(y)), y*)/dnxdny W(y) J(y) h^2	//
			//	K(x(d(y)), y*) = -ik/4 * (term1 + term2 + term3)	//
			//	term1 = nx dot ny * H1(k|x(d(y)) - y*|) / |x(d(y)) - y*|		//
			//	term2 = H0(k|x(d(y))-y*|) * k [ (x-y) dot nx] * [(x-y) dot ny] / |x-y|^2	//
			//	term3 = -2H1(k|x-y|) * [(x-y) dot nx] * [(x-y) dot ny] / |x-y|^3		//
			//	term1R = k/4 { Y1(k|x-y*|)/|x-y*|} * [ nx dot ny ]	//
			//	term1I = -k/4 { J1(k|x-y*|)} * [nx dot ny] / |x-y*|	//
			//	term2R = k/4 {Y0(k|x-y*|) } * k [ (x-y*) dot nx] * [ (x-y*) dot ny] / |x-y*|^2	//
			//	term2I = -k/4 {J0(k|x-y*|) } * k [ (x-y*) dot nx] * [ (x-y*) dot ny]/ |x-y*|^2	//
			//	term3R = k/4 { -2Y1(k|x-y*|) } * [(x-y*) dot nx] * [(x-y*) dot ny] / |x-y*|^3	//
			//	term3I = -k/4 {-2J1(k|x-y*|) } * [(x-y*) dot nx] * [(x-y*) dot ny] / |x-y*|^3	//

			//	x = x(d(y)) = x* + |d(y)| nx	, only outside points	//
			x1 = starxi - fabs(dist[j]) * gradx[i];
			x2 = staryi - fabs(dist[j]) * grady[i];
//			x1 = starxi + fabs(dist[j]) * gradx[i];
//			x2 = staryi + fabs(dist[j]) * grady[i];

			xminusy1 = x1 - starxj;
			xminusy2 = x2 - staryj;
			xminusynorm = sqrt( xminusy1 * xminusy1 + xminusy2 * xminusy2 );
			recipxminusynorm = 1./xminusynorm;

			parameter = wavenum * xminusynorm;

			inner_pd1 = gradx[i] * gradx[j] + grady[i] * grady[j];
			inner_pd2 = xminusy1 * gradx[i] + xminusy2 * grady[i];
			inner_pd3 = xminusy1 * gradx[j] + xminusy2 * grady[j];

			commonterm1 = inner_pd1 * recipxminusynorm;
			commonterm2 = inner_pd2 * recipxminusynorm * inner_pd3 * recipxminusynorm;

			BesselY1value = BesselY1(parameter);
			BesselJ1value = BesselJ1(parameter);

			term1R = BesselY1value * commonterm1;
			term1I = BesselJ1value * commonterm1;
			term2R = BesselY0(parameter) * wavenum * commonterm2;
			term2I = BesselJ0(parameter) * wavenum * commonterm2;
			term3R = -2. * BesselY1value * commonterm2 * recipxminusynorm;
			term3I = -2. * BesselJ1value * commonterm2 * recipxminusynorm;

			if (POLY_TEST)
			{
				Polydoublereal = 0.25 * PolyResult[j] * wavenum * (term1R + term2R + term3R);
				Polydoubleimg = -0.25 * PolyResult[j] * wavenum * (term1I + term2I + term3I);
				PolyKernelH[i][j] = Polysinglereal + Polydoublereal;
				PolyKernelHimg[i][j] = Polysingleimg + Polydoubleimg;
			}
			if (SINE_TEST)
			{
				Sinedoublereal = 0.25 * SineResult[j] * wavenum * (term1R + term2R + term3R);
				Sinedoubleimg = -0.25 * SineResult[j] * wavenum * (term1I + term2I + term3I);
				SineKernelH[i][j] = Sinesinglereal + Sinedoublereal;
				SineKernelHimg[i][j] = Sinesingleimg + Sinedoubleimg;

//				if ( (SineKernelH[i][j] != SineKernelH[i][j]) ||\
//					 (SineKernelHimg[i][j] != SineKernelHimg[i][j]) )
//				{
//					printf("i, j = %d, %d, x(%lf, %lf), y(%lf, %lf)\n", \
//					i, j, x1, x2, starxj, staryj);
//					printf("termsR %lf, %lf, %lf, termsI %lf, %lf, %lf\n", \
//					term1R, term2R, term3R, term1I, term2I, term3I);
//					exit(0);
//				}


			}

//			PolyKernelH[i][j] = Polysinglereal;
//			PolyKernelHimg[i][j] = Polysingleimg;
//			SineKernelH[i][j] = Sinesinglereal;
//			SineKernelHimg[i][j] = Sinesingleimg;

		}
	}
	return 0;
}
