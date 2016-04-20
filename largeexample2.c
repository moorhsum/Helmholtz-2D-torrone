#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//	#include <gsl_sf_bessel.h>
#include <omp.h>
#include "frame.h"
#include "constants.h"
#include "HelmholtzKernel.h"
#include "tools.h"
#include "functions.h"
#include "PDE.h"
#include "BICGSTAB.h"
#include "fd_der.h"
#include "interpolate.h"

//	#define WAVE 5.0	//	determine the wavenumber for (laplace + k^2) u = 0	//

#define ABS_FDTEST 0	//	1 if the finite difference is set to be absolute distance, else a constant multiple of h

#define J_M 1		//	the order of mode for J(z), the Jacobian	//
#define FACTOR 3.0	//	the factor for epsilon = factor * h, this is for the size of the tubular area	//

#define FORCING 1	//	force to print the interior as irrelevant		//
#define DRAW 0		//	print apru (values on the whole frame) or not		//
#define TEST_BESSEL 0	//	test Bessel function for accuracy			//
#define FASTSIZEBOUND 235000	//	for matrix size smaller than this, store the whole thing for speed	//
#define TEST_CONDITION 0	//	print out the matrix for condition number	//
#define TEST_NEW_ESTIMATOR 1	//	test new estimator for solving PDE	//

#define ABS_MODE 0		//	0 for k|nabla d|_1 H, 1 for absolute, 2 for kH	//
#define ABS_EPSILON 0.145
#define ABS_EPSILON0 0.15
#define ABS_DELTA0 0.005
#define COMBO 1
#define DEFAULT_ETA -1.0		//	u = ( dG/dn - i eta G ) * beta , G = (-i/4) H0(k|x|)		//

#define COEFFICIENTA {-7.0, 15.0, 19.0, -14.0}
#define COEFFICIENTB {2.0, 13.0, 16.0, -9.0}

#define OMP_DEFAULT_CACHESIZE 16
#define DEFAULT_OMP_THRESHOLD 1000
#define OMP_DEFAULT_THREADS 12


int TOTALPOINTS = 0;
int TOTALNEWPOINTS = 0;
int N = 256;
int SIZE = 257;
int N2 = 301;
int SIZE2 = 302;
double H, HSQUARED;
double H2;
double WAVE = 15.0;

double NEARFACTOR = 2.0, FARFACTOR = 3.0;
double EPS_ORDER = 0.5, EPSFACTOR = 6.5, DELTAFACTOR = 0.5, DELTA_ORDER = 1.0;
double EPSFACTOR2 = 1.0, DELTAFACTOR2 = 0.0625;
double ABS_NEARDIST = 0.05, ABS_FARDIST = 0.1;

double *Density, *Densityimg, *newPolyDensity, *newPolyDensityimg, *newSineDensity, *newSineDensityimg;
double **newPolyKernelH, **newPolyKernelHimg, **newSineKernelH, **newSineKernelHimg;
double *newzx, *newzy, *newzstarx, *newzstary, *newdist, *newgradx, *newgrady, *newcurv, *newJacobian, *newregfactor;
double *PolyW, *SineW, *PolyResult, *SineResult;
double *newBoundaryValue, *newBoundaryValueimg;

double **DIST, **DIST_DX, **DIST_DY, **CURV;


double get_C(int size, int pi, int pj);
double get_Ccomp( int size, int pi, int pj, double epsilon, double tempC[]);
double omp_get_Ccomp(int size, int pi, int pj, double epsilon, double tempC[]);
int BiCGSTAB(int size);
int BiCGcomplex(int size, double epsilon);
int parallel_BiCGcomplex(int size, double epsilon, int chunk);
int BiCGSTABcomplex(int size, double epsilon);
//	int print_Creg(int size);

int print_f(int size, char *filename, char *filenameimg, double delta0, double epsilon0, int mode );
double cal_u(double x, double y);
double cal_partialu(double x, double y);
double cal_uimg(double x, double y);
double cal_partialuimg(double x, double y);
double cal_vdl(double x, double y);
double cal_vdlimg(double x, double y);

int print_kernel(int size, double samplex, double sampley);
int print_SINGLE_kernel(int size, double samplex, double sampley);
int print_apru(int size);
int NeumannTest(int size, double delta0, double epsilon0, double samplex, double sampley);
int New_NeumannTest(int size, double samplex, double sampley);
int Combo_NeumannTest(int size, double samplex, double sampley, double eta);


int main(int argc, char *argv[])
{
	int i, j;
	int RESOLVE = 0, RESOLVE_NEW = 0;
	int adjust;
	int chunk, nthreads;
	double delta0, epsilon0, eta;
	double thisdist, thisgradx, thisgrady;

	FILE *fpdensity, *fpdensityimg, *fpnewPolydensity, *fpnewPolydensityimg, *fpnewSinedensity, *fpnewSinedensityimg;
	FILE *fpPolyResult, *fpSineResult, *fpnewpointdata;
	FILE *fptestPoly, *fptestSine, *fptestx, *fptestb, *fptestPolyimg, *fptestSineimg, *fptestximg, *fptestbimg;
	FILE *fpdist;

	int counter = 0, newcounter = 0;
	double zx, zy;
//	double h = (INT_U - INT_L)/N;
	double inJ0, outJ0, inY0, outY0, inJ1, outJ1, inY1, outY1;
	double epsilon = FACTOR * H;
	double samplex, sampley;
	char distfilename[50];
	char denfilename[50], denimgfilename[50], newPolydenfilename[50], newPolydenimgfilename[50];
	char newSinedenfilename[50], newSinedenimgfilename[50], PolyResultfilename[50], SineResultfilename[50];
	char newpointdatafilename[50];
	char filenameb[50], filenamebimg[50];
	time_t begin, end, seconds, new_begin, new_end;

	H = (double)(INT_U-INT_L)/(double)(N);
	HSQUARED = H*H;
	H2 = (double)(INT_U-INT_L)/(double)(N2);

	eta = DEFAULT_ETA;

	chunk = OMP_DEFAULT_CACHESIZE;
	nthreads = OMP_DEFAULT_THREADS;
	omp_set_num_threads(nthreads);

	#pragma omp parallel
	{
		#pragma omp master
		{
			printf("Totally using %d threads\n", omp_get_num_threads());
		}
	}



	if (argc > 1)	
	{
		N = atoi(argv[1]);
		SIZE = (N+1);
		H = (INT_U-INT_L)/(double)N;
		HSQUARED = H*H;
		if ( LAPLACE == 0 )
		{
			if (argc > 2)
			{
				WAVE = atof(argv[2]);
			}
			else
			{		
				printf("Enter Wave number k for the helmoholtz problem.\n");
				scanf("%lf", &WAVE);
			}
		}
	}




	if ( (LAPLACE == 0) && (N*PI/(2.0*WAVE) < 6.0) )
	{
		printf("Capturing one wave with roughly %lf points, which may not be enough and result in large error.\n", N*PI/(2.0*WAVE));
		printf("Ideally each wave should have more than six points. In this case, we need at least N = %lf.\n", 12.0*WAVE/PI);
		printf("Adjusting N?(1 for yes, 0 for no) : ");
		scanf("%d", &adjust);

		while (adjust!=0)
		{
			if (adjust==1)
			{
				printf("Please enter new N.\n");
				scanf("%d", &N);
				printf("Using N = %d.\n", (int)(N));
				if (N*PI/WAVE < 12.0)
				{
					printf("The new N = %d still may not be enough.\n", (int)(N));
					printf("Ideally each wave should have more than six points. \n");
					printf("In this case, we need at least N = %lf.\n", 12.0*WAVE/PI);
					printf("Adjusting N? (1 for yes/0 for no) : ");
					scanf("%d", &adjust);
					continue;
				}
				break;
			}
			else
			{
				printf("Invalid input. 1/0 only. Current N = %d\n", (int) (N));
				printf("Ideally each wave should have more than six points. \n");
				printf("In this case, we need at least N = %lf.\n", 12.0*WAVE/PI);
				printf("Adjusting N?(1 for yes/0 for no) : ");
				scanf("%d", &adjust);
			}
		}
	}

	if (TEST_BESSEL==1)
	{
		inJ0 = 0.0;
		printf("Input a number to test BesselJ0 functions. negative to exit.\n");
		while (inJ0 >= -0.0000001)
		{
			scanf("%lf", &inJ0);
			outJ0 = BesselJ0(inJ0);
			printf("Bessel function J0(%lf) is %.10lf.\n", inJ0, outJ0);
		}
	
		inJ1 = 0.0;
		printf("Input a number to test BesselJ1 functions. negative to exit.\n");
		while (inJ1 >= -0.0000001)
		{
			scanf("%lf", &inJ1);
			outJ1 = BesselJ1(inJ1);
			printf("Bessel function J1(%lf) is %.10lf.\n", inJ1, outJ1);
		}

		inY0 = 0.0;
		printf("Input a number to test BesselY0 functions. negative to exit.\n");
		while (inY0 >= -0.0000001)
		{
			scanf("%lf", &inY0);
			outY0 = BesselY0(inY0);
			printf("Bessel function Y0(%lf) is %.10lf.\n", inY0, outY0);
		}

		inY1 = 0.0;
		printf("Input a number to test BesselY1 functions. negative to exit.\n");
		while (inY1 >= -0.0000001)
		{
			scanf("%lf", &inY1);
			outY1 = BesselY1(inY1);
			printf("Bessel function Y1(%lf) is %.10lf.\n", inY1, outY1);
		}
		return 0;
	}

	if ( ( (DIRICHLET == 0) && (LAPLACE == 1) ) || ( (LAPLACE == 0) && (INTERIOR_MODE == 1) ) )
		printf("Check mode. Dirichlet = %d, Laplce = %d, Interior mode = %d.\n", DIRICHLET, LAPLACE, INTERIOR_MODE);

	if (argc < 5)
	{
		printf("Please enter the sample point to solve.\n");
		scanf("%lf%lf", &samplex, &sampley);
	}
	else
	{
		samplex = atof(argv[3]);
		sampley = atof(argv[4]);
		if (argc == 7)
		{
			if (ABS_FDTEST == 1)
			{
				ABS_NEARDIST = atof(argv[5]);
				ABS_FARDIST = atof(argv[6]);
			}
			else
			{
				NEARFACTOR = atof(argv[5]);
				FARFACTOR = atof(argv[6]);
			}
		}
	}

	if (ABS_MODE == 1)
	{
		epsilon = ABS_EPSILON;
		epsilon0 = ABS_EPSILON0;
		delta0 = ABS_DELTA0;
	}
	else
	{
		epsilon = FACTOR * H;
		if (fabs(EPS_ORDER - 1.) < 1e-6)
			epsilon0 = EPSFACTOR * H;
		else
			epsilon0 = EPSFACTOR2 * pow(H, EPS_ORDER);
		if (fabs(DELTA_ORDER - 1.) < 1e-6 )
			delta0 = DELTAFACTOR * H;
		else
			delta0 = DELTAFACTOR2 * pow(H, DELTA_ORDER);
	}

	printf("Using N = %d.  ", (int) (N));
	if ( LAPLACE == 0 )
	{
		if (INTERIOR_MODE==1)
			printf("Solving interior Helmholtz equation where k = %lf.\n", WAVE);
		else
			printf("Solving exterior Helmholtz equation where k = %lf.\n", WAVE);
		printf("This means using roughly %d points to capture a wave.\n", (int)(N*PI/(2.0*WAVE)));
	}
	else
	{
		if (INTERIOR_MODE==1)
			printf("Solving interior Laplace equation.\n");
		else
			printf("Solving exterior Laplace equation.\n");
	}

	counter = 0;
	newcounter = 0;

	if (WAVE < 8.0)
		eta = -2./sqrt( PI*PI + 4.* pow(log(WAVE/2.) + EULER, 2.0));
	else
		eta = -1. * WAVE;

	printf("Eta is equal to %lf. ", eta);

	if (ABS_MODE == 1)
	{
//		epsilon = ABS_EPSILON;
		delta0 = ABS_DELTA0;
		epsilon0 = ABS_EPSILON0;
	}
	else
	{
		if ( fabs(EPS_ORDER - 1.) < 1e-6 )
		{
//			epsilon = FACTOR * H;
			epsilon0 = EPSFACTOR * H;
		}
		else
		{
//			epsilon = FACTOR2 * pow(H, EPS_ORDER);
			epsilon0 = EPSFACTOR2 * pow(H, EPS_ORDER);
		}
		if ( fabs(DELTA_ORDER - 1.) < 1e-6)
		{
			delta0 = DELTAFACTOR * H;
		}
		else
		{
			delta0 = DELTAFACTOR2 * pow(H, DELTA_ORDER);
		}
	}

	if ( (TRAD_ONESIDED == 0) && (COMBO_ONESIDED == 1) )
		epsilon = (epsilon0 - delta0)/2.0;
	else
		epsilon = epsilon0 - delta0;

	if (TEST_NEW_ESTIMATOR == 1)
	{
		printf("Delta = %lfH^%.2lf, epsilon = %lfH^%.2lf.\n", DELTAFACTOR, DELTA_ORDER, EPSFACTOR, EPS_ORDER);
		printf("H = %lf, [delta0, epsilon0] = [%lf, %lf], epsilon = %lf \n", H, delta0, epsilon0, epsilon);
	}

	DIST = (double **) malloc(SIZE * sizeof(double *));
	DIST_DX = (double **) malloc(SIZE * sizeof(double *));
	DIST_DY = (double **) malloc(SIZE * sizeof(double *));
	CURV = (double **) malloc(SIZE * sizeof(double *));

	for (i = 0; i < SIZE; i++)
	{
		DIST[i] = (double *) malloc(SIZE * sizeof(double));
		DIST_DX[i] = (double *) malloc(SIZE * sizeof(double));
		DIST_DY[i] = (double *) malloc(SIZE * sizeof(double));
		CURV[i] = (double *) malloc(SIZE * sizeof(double));
	}
	


	if (DIST_FILEMODE == 1)
	{
		sprintf(distfilename, "distance%d.txt", N);
		if (( fpdist = fopen(distfilename, "r")) == NULL)
		{
			printf("Open %s error.\n", distfilename);
			exit(0);
		}
		for (i = 0; i < SIZE; i++)
		{
			for (j = 0; j < SIZE; j++)
			{
				fscanf(fpdist, "%lf", &DIST[i][j]);
			}
		}

	}
	else
	{
		for (i = 0; i < SIZE; i++)
		{
			zy = INT_L + (double)(i)*H;
			for (j =0 ; j < SIZE; j++)
			{
				zx = INT_L + (double)(j) * H;
				DIST[i][j] = distance(zx, zy);
			}
		}
	}

	#pragma omp parallel for private(i,j) schedule(static, chunk)
	for (i = 0; i < SIZE; i++)
	{
		double Dxp, Dxn, Dyp, Dyn, DDx, DDy;
		double aplus, aminus, bplus, bminus, cplus, cminus, dplus, dminus;
		for (j = 0; j < SIZE; j++)
		{
			FD_DX(DIST, i, j, &Dxp, &Dxn, DER_MODE);
			FD_DY(DIST, i, j, &Dyp, &Dyn, DER_MODE);
			aplus = max(Dxp, 0.0);
			aminus = min(Dxp, 0.0);
			bplus = max(Dxn, 0.0);
			bminus = min(Dxn, 0.0);
			cplus = max(Dyp, 0.0);
			cminus = min(Dyp, 0.0);
			dplus = max(Dyn, 0.0);
			dminus = min(Dyn, 0.0);
			if (DIST[i][j] > 0)
			{
				DIST_DX[i][j] = (fabs(aminus)>fabs(bplus))?aminus:bplus;
				DIST_DY[i][j] = (fabs(cminus)>fabs(dplus))?cminus:dplus;
			}
			else
			{
				DIST_DX[i][j] = (fabs(bminus)>fabs(aplus))?bminus:aplus;
				DIST_DY[i][j] = (fabs(dminus)>fabs(cplus))?dminus:cplus;
			}
			DDx = XDER2(DIST, i, j, 2);
			DDy = XDER2(DIST, i, j, 2);
			CURV[i][j] = -1. * (DDx + DDy);
		}
	}
	for (i = 0; i < SIZE; i++)
	{
		for (j = 0; j < SIZE; j++)
		{
			if ( (DIST[i][j] != DIST[i][j]) || (DIST_DX[i][j] != DIST_DX[i][j]) || \
				(DIST_DY[i][j] != DIST_DY[i][j]) || (CURV[i][j]!= CURV[i][j]) )
			{
				printf("%d, %d dist %lf, dx %lf, dy %lf, curv %lf\n", \
				i, j, DIST[i][j], DIST_DX[i][j], DIST_DY[i][j], CURV[i][j]);
			}
		}
	}



	for ( i = 0; i < SIZE; i++)
	{
		zy = INT_L + H * (double)(i);
		for ( j = 0;j < SIZE; j++)
		{
			zx = INT_L + H * (double)(j);

//			if (DIST_FILEMODE == 1)
//				thisdist = DIST[i][j];
//			else
//				thisdist = distance(zx, zy);

			thisdist = DIST[i][j];


			if (TRAD_ONESIDED == 0)
			{
				if (fabs(thisdist) <= epsilon )
					counter++;
			}
			else
			{			
				if ( ( thisdist >= 0.) && (thisdist <= epsilon) )
				{
//					if (counter > 1200)
//						printf("zx = %.10lf, zy = %.10lf.\n", zx, zy);
					counter++;
				}
			}

			if (COMBO_ONESIDED == 0)
			{
				if ( ( fabs(thisdist) >= delta0) && ( fabs(thisdist) <= epsilon0 ) )
					newcounter++;
			}
			else
			{
				if ( (thisdist >= delta0) && ( thisdist <= epsilon0 ) )
				{
//				if (newcounter > 1200)
//					printf("zx = %.10lf, zy = %.10lf.\n", zx, zy);
					newcounter++;
				}
			}
		}
	}
	printf("The matrix size is %d by %d, new is %d by %d.\n", counter, counter, newcounter, newcounter);
	TOTALPOINTS = counter;
	TOTALNEWPOINTS = newcounter;

	begin = time(NULL);

	sprintf(filenameb, "largerealb.txt");
	sprintf(filenamebimg, "largeimgb.txt");
	print_f(counter, filenameb, filenamebimg, 0., epsilon, 0);

//	printf("Here.\n");

	newBoundaryValue = (double *) malloc(newcounter * sizeof(double));
	newBoundaryValueimg = (double *) malloc(newcounter * sizeof(double));
	if (TEST_NEW_ESTIMATOR == 1)
	{
		sprintf(filenameb, "newlargerealb.txt");
		sprintf(filenamebimg, "newlargeimgb.txt");
		print_f(newcounter, filenameb, filenamebimg, delta0, epsilon0, 1);

	}

	end = time(NULL);
	seconds = end - begin;
//	printf("Printing Boundary condition took %ld minutes %ld seconds.\n", seconds/60, seconds%60);


	Density = (double *) malloc(counter * sizeof(double));
	Densityimg = (double *) malloc(counter * sizeof(double));

	RESOLVE = 0;

	if (LAPLACE == 1)
	{
		sprintf(denfilename, "DensityN%d.txt", (int)(N));
		sprintf(newPolydenfilename, "newPolyDensityN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
		sprintf(newSinedenfilename, "newSineDensityN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
		sprintf(PolyResultfilename, "PolyResultN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
		sprintf(SineResultfilename, "SineResultN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
		sprintf(newpointdatafilename, "newPointdataN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
	}
	else
	{
		sprintf(denfilename, "DensityN%dk%d.txt", (int)(N), (int)(WAVE));
		sprintf(newPolydenfilename, "newPolyDensityN%dk%ddel%.1lfeps%.1lf.txt", \
			(int)(N), (int)(WAVE), DELTAFACTOR, EPSFACTOR);
		sprintf(newSinedenfilename, "newDensityN%dk%ddel%.1lfeps%.1lf.txt", \
			(int)(N), (int)(WAVE), DELTAFACTOR, EPSFACTOR);
		sprintf(PolyResultfilename, "PolyResultN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
		sprintf(SineResultfilename, "SineResultN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
		sprintf(newpointdatafilename, "newPointdataN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
	}


	if (TEST_CONDITION == 1)
	{
//		BiCGcomplex(counter);
		printf("Print matrix for condition number.\n");
		newzx = (double *) malloc(newcounter * sizeof(double));
		newzy = (double *) malloc(newcounter * sizeof(double));
		newzstarx = (double *) malloc(newcounter * sizeof(double));
		newzstary = (double *) malloc(newcounter * sizeof(double));
		newdist = (double * ) malloc(newcounter * sizeof(double));
		newgradx = (double *) malloc(newcounter * sizeof(double));
		newgrady = (double *) malloc(newcounter * sizeof(double));
		newcurv = (double *) malloc(newcounter * sizeof(double));
		newJacobian = (double *) malloc(newcounter * sizeof(double));
		newregfactor = (double *) malloc(newcounter * sizeof(double));
		if (POLY_TEST == 1)
		{
			newPolyDensity = (double *) malloc(newcounter * sizeof(double));
			newPolyDensityimg = (double *) malloc(newcounter * sizeof(double));
			PolyW = (double *) malloc(newcounter * sizeof(double));
			PolyResult = (double *) malloc(newcounter * sizeof(double));
			newPolyKernelH = (double **) malloc(newcounter * sizeof(double *));
			newPolyKernelHimg = (double **) malloc(newcounter * sizeof(double *));
			for ( i = 0; i < newcounter; i++)
			{
				newPolyKernelH[i] = (double *) malloc(newcounter * sizeof(double));
				newPolyKernelHimg[i] = (double *) malloc(newcounter * sizeof(double));
			}
		}
		if (SINE_TEST == 1)
		{
			newSineDensity = (double *) malloc(newcounter * sizeof(double));
			newSineDensityimg = (double *) malloc(newcounter * sizeof(double));
			SineW = (double *) malloc(newcounter * sizeof(double));
			SineResult = (double *) malloc(newcounter * sizeof(double));
			newSineKernelH = (double **) malloc(newcounter * sizeof(double *));
			newSineKernelHimg = (double **) malloc(newcounter * sizeof(double *));
			for ( i = 0; i < newcounter; i++)
			{
				newSineKernelH[i] = (double *) malloc(newcounter * sizeof(double));
				newSineKernelHimg[i] = (double *) malloc(newcounter * sizeof(double));
			}
		}


		if (COMBO == 0)
		{
			new_HelmholtzKernel(newcounter, delta0, epsilon0, newzx, newzy, newzstarx, newzstary, \
			newdist, newgradx, newgrady, newcurv, newJacobian, PolyW, SineW, PolyResult, SineResult, \
			newPolyKernelH, newPolyKernelHimg, newSineKernelH, newSineKernelHimg);
		}
		else
		{
			new_HelmholtzKernelCombo(newcounter, delta0, epsilon0, DIST, DIST_DX, DIST_DY, CURV, \
			newzx, newzy, newzstarx, newzstary, \
			newdist, newgradx, newgrady, newcurv, newJacobian, newregfactor, \
			PolyW, SineW, PolyResult, SineResult, \
			newPolyKernelH, newPolyKernelHimg, newSineKernelH, newSineKernelHimg, WAVE, eta);
		}

		fptestb = fopen("b.txt", "w");
		fptestbimg = fopen("bimg.txt", "w");
		for ( i = 0; i < newcounter; i++)
		{
			fprintf(fptestb, "%.12lf\n", newBoundaryValue[i]);
			fprintf(fptestbimg, "%.12lf\n", newBoundaryValueimg[i]);
		}
		fclose(fptestb);
		fclose(fptestbimg);

		if (POLY_TEST == 1)
		{
			fptestPoly = fopen("Poly.txt", "w");
			fptestPolyimg = fopen("Polyimg.txt", "w");
			for (i = 0; i < newcounter; i++)
			{
				for (j = 0; j < newcounter; j++)
				{
					fprintf(fptestPoly, "%lf\t", newPolyKernelH[i][j]);
					fprintf(fptestPolyimg, "%lf\t", newPolyKernelHimg[i][j]);
				}
				fprintf(fptestPoly, "\n");
				fprintf(fptestPolyimg, "\n");
			}
			free(newPolyDensity);
			free(newPolyDensityimg);
			free(PolyW);
			free(PolyResult);
			for ( i = 0; i < newcounter; i++)
			{
				free(newPolyKernelH[i]);
				free(newPolyKernelHimg[i]);
			}
			free(newPolyKernelH);
			free(newPolyKernelHimg);
			fclose(fptestPoly);
			fclose(fptestPolyimg);
		}
		if (SINE_TEST == 1)
		{
			fptestSine = fopen("Sine.txt", "w");
			fptestSineimg = fopen("Sineimg.txt", "w");

			for ( i = 0; i < newcounter; i++)
			{
				for ( j = 0; j < newcounter; j++)
				{
					fprintf(fptestSine, "%lf\t", newSineKernelH[i][j]);
					fprintf(fptestSineimg, "%lf\t", newSineKernelHimg[i][j]);
				}
				fprintf(fptestSine, "\n");
				fprintf(fptestSineimg, "\n");
			}
			free(newSineDensity);
			free(newSineDensityimg);
			free(SineW);
			free(SineResult);
			for ( i = 0; i < newcounter; i++)
			{
				free(newSineKernelH[i]);
				free(newSineKernelHimg[i]);
			}
			free(newSineKernelH);
			free(newSineKernelHimg);
			fclose(fptestSine);
			fclose(fptestSineimg);
		}
		free(newzx);
		free(newzy);
		free(newzstarx);
		free(newzstary);
		free(newdist);
		free(newgradx);
		free(newgrady);
		free(newcurv);
		free(newJacobian);
		free(newregfactor);
		BiCGcomplex(counter, epsilon);
	}

//////////////////////////////////////////////////////////////	BEGIN OF SINGLE LAYER APPROACH	///////////////

	if (TEST_SINGLE_LAYER == 1)
	{

	printf("\nSingle Layer Approach:\t");

	if ((fpdensity = fopen( denfilename, "r"))==NULL)
	{
//		printf("Open file %s error.\n", denfilename);
//		if (LAPLACE == 1)
//			printf("File %s does not exist.\n", denfilename);
		RESOLVE = 1;
	}
	if (LAPLACE == 0)
	{
		sprintf(denimgfilename, "DensityimgN%dk%d.txt", (int)(N), (int)(WAVE));
		if ((fpdensityimg = fopen(denimgfilename, "r"))==NULL)
		{
//			printf("Open file %s error.\n", denimgfilename);
//			printf("File %s does not exist.\n", denimgfilename);
			RESOLVE = 1;
		}
	}
	if (RESOLVE == 1)
	{
		printf("Resolving double layer density...\n");
		begin = time(NULL);
		if (LAPLACE == 0)
		{
			if ( (counter < DEFAULT_OMP_THRESHOLD) || (counter > FASTSIZEBOUND) )
				BiCGcomplex(counter, epsilon);
			else
				parallel_BiCGcomplex(counter, epsilon, chunk);
		}
//		BiCGSTABcomplex(counter);
		else
			BiCGSTAB(counter);

		end = time(NULL);
		seconds = end - begin;
		printf("BiCG took %ld minutes %ld seconds.\n", seconds/60, seconds%60);
	}
	else
	{
		printf("Obtaining double layer density from file...\n");
		for ( i = 0; i < counter; i++)
		{
			fscanf(fpdensity, "%lf", &Density[i]);
			if (LAPLACE == 0)
				fscanf(fpdensityimg, "%lf", &Densityimg[i]);
		}
		fclose(fpdensity);
		if (LAPLACE == 0)
			fclose(fpdensityimg);
	}


	begin = time(NULL);
	if (PRINT_KERNEL == 1)
	{
		if (DIRICHLET == 1)
			print_kernel(counter, samplex, sampley);
		else
			print_SINGLE_kernel(counter, samplex, sampley);
	}
	if (DRAW==1)
		print_apru(counter);
//	end = time(NULL);
//	seconds = end-begin;
//	printf("Printing sample kernel took %ld minutes %ld seconds.\n", seconds/60, seconds%60);
	
	NeumannTest(counter, 0., epsilon, samplex, sampley);

	free(Density);
	free(Densityimg);

	}

/////////////////////////////////////////	END OF SINGLE LAYER APPROACH, BEGIN OF COMBINATION APPROACH	/////////

	if (TEST_NEW_ESTIMATOR == 1)
	{
		printf("\nCombination Approach:\t");
		if ( (fpnewpointdata = fopen( newpointdatafilename, "r")) == NULL)
		{
//			printf("File %s does not exist.\n", newpointdatafilename);
			fpnewpointdata = fopen(newpointdatafilename, "w");
			RESOLVE_NEW = 1;
		}
		if (POLY_TEST == 1)
		{
			if ( (fpnewPolydensity = fopen( newPolydenfilename, "r"))==NULL)
			{
//				if (LAPLACE == 1)
//					printf("File %s does not exist.\n", newPolydenfilename);
				fpnewPolydensity = fopen( newPolydenfilename, "w");
				RESOLVE_NEW = 1;
			}
			if (LAPLACE == 0)
			{
				sprintf(newPolydenimgfilename, "newPolyDensityimgN%dk%ddel%.1lfeps%.1lf.txt", \
				(int)(N), (int)(WAVE), DELTAFACTOR, EPSFACTOR);
				if ( (fpnewPolydensityimg = fopen(newPolydenimgfilename, "r")) == NULL)
				{
//					printf("File %s does not exist.\n", newPolydenimgfilename);
					fpnewPolydensityimg = fopen(newPolydenimgfilename, "w");
					RESOLVE_NEW = 1;
				}
			}
			if ( (fpPolyResult = fopen( PolyResultfilename, "r")) == NULL)
			{
//				printf("File %s does not exist.\n", PolyResultfilename);
				fpPolyResult = fopen( PolyResultfilename, "w");
				RESOLVE_NEW = 1;
			}
		}
		if (SINE_TEST == 1)
		{
			if ( (fpnewSinedensity = fopen( newSinedenfilename, "r"))==NULL)
			{
//				if (LAPLACE == 1)
//					printf("File %s does not exist.\n", newSinedenfilename);
				fpnewSinedensity = fopen( newSinedenfilename, "w");
				RESOLVE_NEW = 1;
			}
			if (LAPLACE == 0)
			{
				sprintf(newSinedenimgfilename, "newSineDensityimgN%dk%ddel%.1lfeps%.1lf.txt", \
				(int)(N), (int)(WAVE), DELTAFACTOR, EPSFACTOR);
				if ( (fpnewSinedensityimg = fopen(newSinedenimgfilename, "r")) == NULL)
				{
//					printf("File %s does not exist.\n", newSinedenimgfilename);
					fpnewSinedensityimg = fopen(newSinedenimgfilename, "w");
					RESOLVE_NEW = 1;
				}
			}
			if ( (fpSineResult = fopen( SineResultfilename, "r")) == NULL)
			{
//				printf("File %s does not exist.\n", SineResultfilename);
				fpSineResult = fopen( SineResultfilename, "w");
				RESOLVE_NEW = 1;
			}
		}	

		newzx = (double *) malloc(newcounter * sizeof(double));
		newzy = (double *) malloc(newcounter * sizeof(double));
		newzstarx = (double *) malloc(newcounter * sizeof(double));
		newzstary = (double *) malloc(newcounter * sizeof(double));
		newdist = (double * ) malloc(newcounter * sizeof(double));
		newgradx = (double *) malloc(newcounter * sizeof(double));
		newgrady = (double *) malloc(newcounter * sizeof(double));
		newcurv = (double *) malloc(newcounter * sizeof(double));
		newJacobian = (double *) malloc(newcounter * sizeof(double));
		newregfactor = (double *) malloc(newcounter * sizeof(double));
		if (POLY_TEST == 1)
		{
			newPolyDensity = (double *) malloc(newcounter * sizeof(double));
			newPolyDensityimg = (double *) malloc(newcounter * sizeof(double));
			PolyW = (double *) malloc(newcounter * sizeof(double));
			PolyResult = (double *) malloc(newcounter * sizeof(double));
		}
		if (SINE_TEST == 1)
		{
			newSineDensity = (double *) malloc(newcounter * sizeof(double));
			newSineDensityimg = (double *) malloc(newcounter * sizeof(double));
			SineW = (double *) malloc(newcounter * sizeof(double));
			SineResult = (double *) malloc(newcounter * sizeof(double));
		}
		if (RESOLVE_NEW == 1)
		{
			printf("Resolving combo density...\n");
			new_begin = time(NULL);
			if (POLY_TEST == 1)
			{
				newPolyKernelH = (double **) malloc(newcounter * sizeof(double *));
				newPolyKernelHimg = (double **) malloc(newcounter * sizeof(double *));
				for ( i = 0; i < newcounter; i++)
				{
					newPolyKernelH[i] = (double *) malloc(newcounter * sizeof(double));
					newPolyKernelHimg[i] = (double *) malloc(newcounter * sizeof(double));
				}
			}
			if (SINE_TEST == 1)
			{
				newSineKernelH = (double **) malloc(newcounter * sizeof(double *));
				newSineKernelHimg = (double **) malloc(newcounter * sizeof(double *));
				for ( i = 0; i < newcounter; i++)
				{
					newSineKernelH[i] = (double *) malloc(newcounter * sizeof(double));
					newSineKernelHimg[i] = (double *) malloc(newcounter * sizeof(double));
				}
			}
			if (LAPLACE == 0)
			{

				if (COMBO == 0)
				{
				new_HelmholtzKernel(newcounter, delta0, epsilon0, newzx, newzy, newzstarx, newzstary, \
				newdist, newgradx, newgrady, newcurv, newJacobian, PolyW, SineW, PolyResult, SineResult, \
				newPolyKernelH, newPolyKernelHimg, newSineKernelH, newSineKernelHimg);
				}
				else
				{
				new_HelmholtzKernelCombo(newcounter, delta0, epsilon0, DIST, DIST_DX, DIST_DY, CURV, \
				newzx, newzy, newzstarx, newzstary,\
				newdist, newgradx, newgrady, newcurv, newJacobian, newregfactor, \
				PolyW, SineW, PolyResult, SineResult, \
				newPolyKernelH, newPolyKernelHimg, newSineKernelH, newSineKernelHimg, WAVE, eta);
				}

////////////////////////////////////////////////////////////////////////////////////////
/*
				fptestA = fopen("A.txt", "w");
				fptestb = fopen("b.txt", "w");
				fptestAimg = fopen("Aimg.txt", "w");
				fptestbimg = fopen("bimg.txt", "w");
				for (int i = 0; i < newcounter; i++)
				{
					for (int j = 0; j < newcounter; j++)
					{
						fprintf(fptestA, "%.12lf\t", newPolyKernelH[i][j]);
						fprintf(fptestAimg, "%.12lf\t", newPolyKernelHimg[i][j]);
					}
					fprintf(fptestb, "%.12lf\n", newBoundaryValue[i]);
					fprintf(fptestbimg, "%.12lf\n", newBoundaryValueimg[i]);
					fprintf(fptestA, "\n");
					fprintf(fptestAimg, "\n");
				}
				fclose(fptestA);
				fclose(fptestAimg);
				fclose(fptestb);
				fclose(fptestbimg);
*/
//////////////////////////////////////////////////////////////////////////////////////////

				for ( i = 0; i < newcounter; i++)
				{
					fprintf(fpnewpointdata, "%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\n", \
					newdist[i], newzstarx[i], newzstary[i], newgradx[i], newgrady[i] );
				}
				fclose(fpnewpointdata);

				if (POLY_TEST == 1)
				{
					if (newcounter < DEFAULT_OMP_THRESHOLD )
					{
					CompBiCG(newcounter, newPolyDensity, newPolyDensityimg, newPolyKernelH, \
					newPolyKernelHimg, newBoundaryValue, newBoundaryValueimg);
					}
					else
					{
					omp_CompBiCG(newcounter, newPolyDensity, newPolyDensityimg, newPolyKernelH, \
					newPolyKernelHimg, newBoundaryValue, newBoundaryValueimg, chunk);
					}
					for ( i = 0; i < newcounter; i++)
					{
						fprintf(fpnewPolydensity, "%.12lf\n", newPolyDensity[i]);
						fprintf(fpnewPolydensityimg, "%.12lf\n", newPolyDensityimg[i]);
						fprintf(fpPolyResult, "%.12lf\n", PolyResult[i]);

						free(newPolyKernelH[i]);
						free(newPolyKernelHimg[i]);
					}
					free(newPolyKernelH);
					free(newPolyKernelHimg);
					fclose(fpnewPolydensity);
					fclose(fpnewPolydensityimg);
					fclose(fpPolyResult);
				}

				if (SINE_TEST == 1)
				{
					if (newcounter < DEFAULT_OMP_THRESHOLD )
					{
					CompBiCG(newcounter, newSineDensity, newSineDensityimg, newSineKernelH, \
					newSineKernelHimg, newBoundaryValue, newBoundaryValueimg);
					}
					else
					{
					omp_CompBiCG(newcounter, newSineDensity, newSineDensityimg, newSineKernelH, \
					newSineKernelHimg, newBoundaryValue, newBoundaryValueimg, chunk);
					}

//				fptestx = fopen("x.txt", "w");
//				fptestximg = fopen("ximg.txt", "w");
					for ( i = 0; i < newcounter; i++)
					{
//					fprintf(fptestx, "%.12lf\n", newPolyDensity[i]);
//					fprintf(fptestximg, "%.12lf\n", newPolyDensityimg[i]);
						fprintf(fpnewSinedensity, "%.12lf\n", newSineDensity[i]);
						fprintf(fpnewSinedensityimg, "%.12lf\n", newSineDensityimg[i]);
						fprintf(fpSineResult, "%.12lf\n", SineResult[i]);

						free(newSineKernelH[i]);
						free(newSineKernelHimg[i]);
					}
					free(newSineKernelH);
					free(newSineKernelHimg);
//				fclose(fptestx);
//				fclose(fptestximg);
					fclose(fpnewSinedensity);
					fclose(fpnewSinedensityimg);
					fclose(fpSineResult);
				}
			}
			else
			{
				printf("Not now.\n");
				exit(0);
//				newBiCGSTAB
			}
			new_end = time(NULL);
			seconds = new_end-new_begin;
			printf("Solving Combo took %ld minutes %ld seconds.\n", seconds/60, seconds%60);
		}
		else
		{
			printf("Obtaining combo density from file...\n");
			for ( i = 0; i < newcounter; i++)
			{
				fscanf(fpnewpointdata, "%lf%lf%lf%lf%lf", \
				&newdist[i], &newzstarx[i], &newzstary[i], &newgradx[i], &newgrady[i]);
				if (POLY_TEST == 1)
				{
					fscanf(fpnewPolydensity, "%lf", &newPolyDensity[i]);
					if (LAPLACE == 0)
						fscanf(fpnewPolydensityimg, "%lf", &newPolyDensityimg[i]);
					fscanf(fpPolyResult, "%lf", &PolyResult[i]);
				}
				if (SINE_TEST == 1)
				{
					fscanf(fpnewSinedensity, "%lf", &newSineDensity[i]);
					if (LAPLACE == 0)
						fscanf(fpnewSinedensityimg, "%lf", &newSineDensityimg[i]);
					fscanf(fpSineResult, "%lf", &SineResult[i]);
				}
			}
			if (POLY_TEST == 1)
			{
				fclose(fpnewPolydensity);
				fclose(fpPolyResult);
				if (LAPLACE == 0)
					fclose(fpnewPolydensityimg);
			}
			if (SINE_TEST == 1)
			{
				fclose(fpnewSinedensity);
				fclose(fpSineResult);
				if (LAPLACE == 0)
					fclose(fpnewSinedensityimg);
			}
			fclose(fpnewpointdata);
//			printf("Here.\n");
		}

		if (COMBO == 0)
			New_NeumannTest(newcounter, samplex, sampley);
		else
			Combo_NeumannTest(newcounter, samplex, sampley, eta);

		if (POLY_TEST == 1)
		{
			free(newPolyDensity);
			free(newPolyDensityimg);
			free(PolyW);
			free(PolyResult);
		}
		if (SINE_TEST == 1)
		{
			free(newSineDensity);
			free(newSineDensityimg);
			free(SineW);
			free(SineResult);
		}
		free(newzx);
		free(newzy);
		free(newzstarx);
		free(newzstary);
		free(newdist);
		free(newgradx);
		free(newgrady);
		free(newcurv);
		free(newJacobian);
		free(newregfactor);
	}

	return 0;
}


int print_apru(int size)
{
	int i, j, k;
	FILE *fpapru, *fpapruimg, *fpx, *fpximg, *fpinwave, *fpinwaveimg;
	double x[size], ximg[size];
	double zx, zy, tempwx, tempwy, tstarx, tstary, partialx, partialy;
	double partialximg, partialyimg;
	double wx[size], wy[size], wstarx[size], wstary[size];
	double J[size], delta[size], gradx[size], grady[size], result[size];
	double thisdist, thisgradx, thisgrady;

	double sum, sumimg, pfunction, temp, tempimg;
	double epsilon = FACTOR * H;

	double nonhomogeneous[2] = {0.0, 0.0};

	char denfilename[50], denimgfilename[50];

	int counter = 0;

	if ((fpapru = fopen( "apru.txt", "w+"))==NULL)
	{
		printf("apru.txt error.\n");
		exit(0);
	}

	if (LAPLACE == 1)
		sprintf(denfilename, "DensityN%d.txt", (int)(N));
	else
		sprintf(denfilename, "DensityN%dk%d.txt", (int)(N), (int)(WAVE));
		
	if ((fpx = fopen( denfilename, "r"))==NULL)
	{
		printf("Open %s error.\n", denfilename);
		exit(0);
	}

	if (LAPLACE==0)
	{
		if ((fpapruimg = fopen( "apruimg.txt", "w+"))==NULL)
		{
			printf("apruimg.txt error.\n");
			exit(0);
		}
		sprintf(denimgfilename, "DensityimgN%dk%d.txt", (int)(N), (int)(WAVE));
		if ((fpximg = fopen(denimgfilename, "r"))==NULL)
		{
			printf("Open %s error.\n", denimgfilename);
			exit(0);
		}
		if ((fpinwave = fopen("inwave.txt", "w+"))==NULL)
		{
			printf("Open inwave.txt error.\n");
			exit(0);
		}
		if ((fpinwaveimg = fopen("inwaveimg.txt", "w+"))==NULL)
		{
			printf("Open inwaveimg.txt error.\n");
			exit(0);
		}
	}


	for ( i = 0; i < size; i++)
	{
		fscanf(fpx, "%lf", &x[i]);
		if (LAPLACE==0)
			fscanf(fpximg, "%lf", &ximg[i]);
		if (i < 0)
		{
			printf("here in apru.\n");
			printf("x[%d] = %lf, ", i, x[i]);
			printf("ximg[%d] = %lf\n", i, ximg[i]);
		}
	//	printf("%lf ", x[i]);	
	}

	epsilon = FACTOR * H;
	for ( i = 0; i < SIZE; i++)
	{
		tempwy = INT_L + H*(double)(i);
		for ( j = 0; j < SIZE; j++)
		{
			tempwx = INT_L + H*(double)(j);
//			printf("distance = %lf, epsilon = %lf.\n", distance(tempwx, tempwy), epsilon);
//			if (DIST_FILEMODE == 1)
//				thisidst = DIST[i][j];
//			else
//				thisdist = distance(tempwx, tempwy);

			thisdist = DIST[i][j];

			thisgradx = gradientdx(tempwx, tempwy);
			thisgrady = gradientdy(tempwx, tempwy);
//			epsilon = ( fabs(thisgradx) + fabs(thisgrady) ) * FACTOR * H;
			if ( fabs(thisdist) > epsilon )
				continue;

//			printf("zx = %lf, zy = %lf.\n", tempwx, tempwy);
//			tstarx = starx(tempwx, tempwy);
			tstarx = tempwx - thisdist * DIST_DX[i][j];
			tstary = tempwy - thisdist * DIST_DY[i][j];
//			printf("counter = %d.\n", counter);
//			tstary = stary(tempwx, tempwy);
//			printf("counter = %d.\n", counter);
			wx[counter] = tempwx;
			wy[counter] = tempwy;
			wstarx[counter] = tstarx;
			wstary[counter] = tstary;


			delta[counter] = cal_delta(epsilon, thisdist);
//			J[counter] = cal_J(J_M, tempwx, tempwy);
			J[counter] = 1. + thisdist * CURV[i][j];

//			gradx[counter] = cal_interpgrad(1, tstarx, tstary);
//			grady[counter] = cal_interpgrad(2, tstarx, tstary);

			cal_interpgrad(tstarx, tstary, DIST_DX, DIST_DY, &gradx[counter], &grady[counter]);

			result[counter] = H * H * delta[counter] * J[counter];

			if (result[counter] != result[counter])
			{
				printf("counter = %d, delta = %lf.\n", counter, delta[counter]);
				exit(0);
			}

			counter++;
		}
	}
	

//	printf("here.\n");
	if (counter!=size)
	{
		printf("Apru size wrong.\n");
		exit(0);
	}
	counter = 0;
	for (i = 0; i < SIZE2; i++)
	{
		zy = INT_L + H2 * (double)(i);
		for ( j = 0; j < SIZE2; j++)
		{
			zx = INT_L + H2 * (double)(j);
			thisdist = DIST[i][j];
			sum = 0.0;
			sumimg = 0.0;
//			counter = 0;
			if (LAPLACE == 0)
			{
				if ( ((zx > 0.0) && (fabs(zy) < radius)) || (thisdist > -0.015) )
				{
					fprintf(fpinwave, "%.8lf ", 0.0);
					fprintf(fpinwaveimg, "%.8lf ", 0.0);
				}
				else
				{
					fprintf(fpinwave, "%.8lf ", cos(WAVE * zx));
					fprintf(fpinwaveimg, "%.8lf ", sin(WAVE * zx));
				}

				if (FORCING == 1)
				{
					if ( thisdist > -0.015)
					{
						fprintf(fpapru, "%.8lf ", sum);
						fprintf(fpapruimg, "%.8lf ", sum);
						continue;
					}
				}
			}

			for ( k = 0; k < size; k++)
			{
				tstarx = wstarx[k];
				tstary = wstary[k];

				partialx = cal_partial(1, zx, zy, tstarx, tstary);
				partialy = cal_partial(2, zx, zy, tstarx, tstary);
				partialximg = cal_partial(3, zx, zy, tstarx, tstary);
				partialyimg = cal_partial(4, zx, zy, tstarx, tstary);
			
			
				//	for interior problem	//	
				if (INTERIOR_MODE==1)
				{
					temp = (-1.0) * (partialx*gradx[k] + partialy*grady[k]); // outer normal
					if (LAPLACE==0)
						tempimg = (-1.0) * (partialximg*gradx[k] + partialyimg*grady[k]);
				}
				//	for exterior problem	//
				else
				{
					temp = partialx*gradx[k] + partialy*grady[k];		// inner normal
					if (LAPLACE==0)
						tempimg = partialximg*gradx[k] + partialyimg*grady[k];
				}
				if (LAPLACE==1)
				{
					sum += (temp * result[k] * x[k]);
				}
				else
				{
					sum += (result[k] * (temp * x[k] - tempimg * ximg[k]));
					sumimg += (result[k] * (tempimg * x[k] + temp * ximg[k]));
				}

				if ((temp!=temp)||(temp+1.0==temp)||(tempimg!=tempimg)||(tempimg+1.0==tempimg))
				{
					printf("zx = %lf, zy = %lf, counter = %d, ", zx, zy, counter);
					printf("wx = %lf, wy = %lf, wstarx = %lf, wstary = %lf.\n", wx[k], wy[k], wstarx[k], wstary[k]);
					printf("delta = %lf, pfunction = %lf, temp = %lf\n", delta[k], pfunction, temp);
				}
//				if (counter < 10)
//					printf("%lf", sum);
//				counter++;

			}

			if (NONHOM_MODE == 1)
			{
				cal_nonhom(zx, zy, nonhomogeneous);
				sum += nonhomogeneous[0];
				sumimg += nonhomogeneous[1];
			}
			fprintf(fpapru, "%.8lf ", sum);
			if (LAPLACE==0)
				fprintf(fpapruimg, "%.8lf ", sumimg);

		}
		fprintf(fpapru, "\n");
		if (LAPLACE == 0)
		{
			fprintf(fpinwave, "\n");
			fprintf(fpinwaveimg, "\n");
			fprintf(fpapruimg, "\n");
		}
	}

	fclose(fpapru);
	fclose(fpx);
	if (LAPLACE==0)
	{
		fclose(fpinwave);
		fclose(fpximg);
		fclose(fpinwaveimg);
		fclose(fpapruimg);
	}
	return 0;
}



int print_f(int size, char *filename, char *filenameimg, double delta0, double epsilon0, int mode)
{
	int i, j;
	double f[size], fimg[size];
	double zx, zy, zstarx, zstary;
//	double epsilon = FACTOR*h;
	double tempvalue[2];
	double thisdist, thisgradx, thisgrady;
	int counter = 0;

	FILE *fpf, *fpfimg;

//	printf("here inside print_f.\n");
	if ((fpf = fopen(filename, "w"))==NULL)
	{
		printf("open file %s error.\n", filename);
		return 0;
	}
	if (LAPLACE == 0)
	{
		if ((fpfimg = fopen(filenameimg, "w"))==NULL)
		{
			printf("Open file %s error.\n", filenameimg);
			return 0;
		}
		fprintf(fpfimg, "%d\n", size);
	}

	fprintf(fpf, "%d\n", size);
	counter = 0;

//	if (ABS_MODE != 1)
//		epsilon = FACTOR * h;

	for ( i = 0; i < SIZE; i++)
	{
		zy = INT_L + H*i;
		for ( j = 0; j < SIZE; j++)
		{
			zx = INT_L + H*j;
//			if (ABS_MODE == 1)
//			{
//				if (mode == 0)
//				{
//					epsilon = ABS_EPSILON;
//					delta = -1.;
//				}
//				else
//				{
//					epsilon = ABS_EPSILON0;
//					delta = ABS_DELTA0;
//				}
//			}
//			else
//			{
//				if ( (mode == 0) && (ABS_MODE == 0) )
//				{
//					thisgradx = gradientdx(zx, zy);
//					thisgrady = gradientdy(zx, zy);
//					epsilon = ( fabs(thisgradx) + fabs(thisgrady) ) * FACTOR * h;
//				}
//			}

//			thisdist = distance(zx, zy);

			thisdist = DIST[i][j];

//			printf("ABS_MODE = %d, counter = %d, dist = %lf, epsilon = %lf, delta = %lf.\n", \
//			ABS_MODE, counter, thisdist, epsilon, delta);


//			if ( (fabs(thisdist) > epsilon) || (fabs(thisdist) < delta) )

			if (mode == 0)
			{
				if (TRAD_ONESIDED == 0)
				{
					if ( (fabs(thisdist) > epsilon0) || (fabs(thisdist) < delta0) )
						continue;
				}
				else
				{
					if ( (thisdist > epsilon0) || (thisdist < delta0) )
						continue;
				}
			}
			else
			{
				if (COMBO_ONESIDED == 0)
				{
					if ( ( fabs(thisdist) > epsilon0) || (fabs(thisdist) < delta0) )
						continue;
				}
				else
				{
					if ( (thisdist > epsilon0) || (thisdist < delta0) )
						continue;
				}
			}

//			zstarx = starx(zx, zy);
//			zstary = stary(zx, zy);

			zstarx = zx - thisdist * DIST_DX[i][j];
			zstary = zy - thisdist * DIST_DY[i][j];

			if (DIRICHLET == 1)
				f[counter] = cal_u(zstarx, zstary);
			else
			{
				f[counter] = cal_partialu(zstarx, zstary);
//				f[counter] = 2. * cal_partialu(zstarx, zstary);
			}

			if (LAPLACE == 0)
			{
				if (DIRICHLET == 1)
					fimg[counter] = cal_uimg(zstarx, zstary);
				else
				{
					fimg[counter] = cal_partialuimg(zstarx, zstary);
//					fimg[counter] = 2. * cal_partialuimg(zstarx, zstary);
				}
			}

			if (NONHOM_MODE==1)
			{
				cal_nonhom(zstarx, zstary, tempvalue);
				f[counter] -= tempvalue[0];	//	integral of psi and phi for inhomogeneous mode	//
				if (LAPLACE == 0)
					fimg[counter] -= tempvalue[1];
			}
			fprintf(fpf, "%.15lf ", f[counter]);
			if (LAPLACE == 0)
				fprintf(fpfimg, "%.15lf ", fimg[counter]);

			if (mode == 1)
			{
				newBoundaryValue[counter] = f[counter];
				newBoundaryValueimg[counter] = fimg[counter];
			}
//			if (counter > 1200)
//			{
//				printf("(i, j) = (%d, %d), counter = %d, dist = %lf, epsilon = %lf, delta = %lf\n", \
//				i, j, counter, thisdist, epsilon, delta);
//			}
//			printf("counter = %d, epsilon = %lf, delta = %lf\n", counter, epsilon, delta);
			counter++;
		}
	}
	fclose(fpf);
	if (LAPLACE == 0)
		fclose(fpfimg);
	return 0;
}


void get_kite_normal(double x, double y, double *Dx, double *Dy)
{
	double sint, cost;

	sint = y/1.5;
	cost = x + 1.3*sint*sint;

	*Dx = -1. * sint - 2.6*sint*cost;
	*Dy = 1.5 * cost;
}


double cal_u(double x, double y)
{
	int i, j;
	const static double a[4] = COEFFICIENTA;
	const static double b[4] = COEFFICIENTB;
	static double Lambda[6][6] = ROOTBESSEL;
	static double Lambda2[6][6] = ROOTBESSEL;
	static int lambda_flag = 0;
	static int Bessel_flag = 0;
	static int kite_flag = 0;
	double r, r2, temptheta, theta;
	double sum = 0.0;
	static double J0prime, Y0prime, Besselnorm;

	double Dx, Dy;
	if (KITE_TEST == 1)
	{
//		get_kite_normal(x,y,&Dx, &Dy);
		if (kite_flag == 0)
			printf("Kite shape.\n");
		kite_flag = 1;
		return 0.0;
//		exit(0);
	}





	if (Bessel_flag == 0)
	{
//		J0prime = bessj0(WAVE*radius);
//		Y0prime = bessy0(WAVE*radius);
		J0prime = BesselJ0(WAVE*radius);
		Y0prime = BesselY0(WAVE*radius);
		
		Besselnorm = J0prime*J0prime + Y0prime*Y0prime;
		Bessel_flag = 1;
	}
	double J0 = 0.0;
	double Y0 = 0.0;

	if ( (LAPLACE==0) && (NONHOM_MODE==1) )		//	Total wave Dirichlet hard condition.	//
		return 0.0;

//	return 1.0;
//	if (LAPLACE==0)
//	{
//		if (x < 0.0)
//			return (-1.0 * cos(WAVE * x) );
//		else
//			return 0.0;
//	}
//	printf("Deep in cal_u.\n");
	
	if (lambda_flag == 0)
	{
		for ( i = 0; i < 6; i++)
			for ( j = 0; j < 6; j++)
			{
				Lambda[i][j] /= radius;
				Lambda2[i][j] /= radius2;
			}
		lambda_flag++;
	}

	r = sqrt(pow(x-Cx,2.0) + pow(y-Cy,2.0));
	r2 = sqrt(pow(x-Cx2,2.0) + pow(y-Cy2,2.0));
	temptheta = atan((y-Cy)/(x-Cx));	// note the range of inverse trig functions is tricky	//
	if ((x >= Cx) && (y < Cy))	//	4th quadrant
		theta = temptheta + 2.0*PI;
	else if ((x >= Cx) && (y >= Cy))	//	1st quadrant
		theta = temptheta;
	else
		theta = temptheta + PI;

////////////////////////////////////////////////////////////////////////

	if (LAPLACE==0)
	{

		if (INTERIOR_MODE == 1)
		{
//			J0 = BesselJ0(WAVE*r);
//			J0prime = BesselJ0(WAVE*radius);
//			Y0 = BesselY0(WAVE*r);
//			Y0prime = BesselY0(WAVE*radius);
//			return ( J0*J0prime + Y0*Y0prime )/Besselnorm;

			return BesselJ0(WAVE*r)/BesselJ0(WAVE*radius);
		}
		else
		{
			return BesselJ0(WAVE*r);

			J0 = BesselJ0(WAVE*r);
//	J0prime = BesselJ0(WAVE*radius);
			Y0 = BesselY0(WAVE*r);
//	Y0prime = BesselY0(WAVE*radius);
			return ( J0*J0prime + Y0*Y0prime )/Besselnorm;
	
//	return 0.0;
//	return BesselJ0(Lambda[0][1]*r);

//	return BesselJ0(Lambda[0][1]*r) + BesselJ0(Lambda2[0][1]*r2);
//	return BesselJ0(Lambda[0][1]*r) + BesselJ0(Lambda2[0][2]*r2) + BesselJ0(Lambda[0][3]*r);
//	return BesselJ0(Lambda[0][3]*r) + 2.0*BesselJ0(Lambda[0][3]*r2);
		}

	}
	else
	{
		if (INTERIOR_MODE==0)
		{
			for ( i = 1; i <= 3; i++)
				sum += ( pow((1.0/r),(double)(i)) * ( a[i-1]*cos(i*theta) + b[i-1]*sin(i*theta) ) );
			return sum;
		}
		else
		{
			for ( i = 1; i <= 4; i++)
				sum += ( pow(r,(double)(i)) * ( a[i-1]*cos(i*theta) + b[i-1]*sin(i*theta) ) );
			return sum;
		}
	}

////////////////////////////////////////////////////////////////////////

	for ( i = 1; i <= 4; i++)
		sum += pow( r, (double)(i) ) * ( a[i-1]*cos(i*theta) + b[i-1]*sin(i*theta) );

	return sum;
}

double cal_uimg(double x, double y)
{
	int i, j;
	const static double a[4] = COEFFICIENTA;
	const static double b[4] = COEFFICIENTB;
	static double Lambda[6][6] = ROOTBESSEL;
	static double Lambda2[6][6] = ROOTBESSEL;
	static int lambda_flag = 0;
	static int Bessel_flag = 0;
	static int kite_flag = 0;
	double r, r2, temptheta, theta;
	double sum = 0.0;
	static double J0prime, Y0prime, Besselnorm;

	double Dx, Dy;
	double J0 = 0.0;
	double Y0 = 0.0;


	if (KITE_TEST == 1)
	{
		if (kite_flag == 0)
			printf("Kite test.\n");
		kite_flag = 1;
		return 0.0;
//		exit(0);
	}

	if ( (LAPLACE==0) && (NONHOM_MODE==1) )		//	Total wave Dirichlet hard condition.	//
		return 0.0;

	if (Bessel_flag == 0)
	{
		J0prime = BesselJ0(WAVE*radius);
		Y0prime = BesselY0(WAVE*radius);
		Besselnorm = J0prime*J0prime + Y0prime*Y0prime;
	}	


//	return 1.0;
//	if (x < 0.0)
//		return (-1.0 * sin(WAVE * x) );
//	else
//		return 0.0;
//	printf("Deep in cal_u.\n");
	
	if (lambda_flag == 0)
	{
		for ( i = 0; i < 6; i++)
			for ( j = 0; j < 6; j++)
			{
				Lambda[i][j] /= radius;
				Lambda2[i][j] /= radius2;
			}
		lambda_flag++;
	}

	r = sqrt(pow(x-Cx,2.0) + pow(y-Cy,2.0));
	r2 = sqrt(pow(x-Cx2,2.0) + pow(y-Cy2,2.0));
	temptheta = atan((y-Cy)/(x-Cx));	// note the range of inverse trig functions is tricky	//
	if ((x >= Cx) && (y < Cy))	//	4th quadrant
		theta = temptheta + 2.0*PI;
	else if ((x >= Cx) && (y >= Cy))	//	1st quadrant
		theta = temptheta;
	else
		theta = temptheta + PI;

////////////////////////////////////////////////////////////////////////

	if (LAPLACE==0)
	{

		if (INTERIOR_MODE == 1)
		{
//			J0 = BesselJ0(WAVE*r);
//			J0prime = BesselJ0(WAVE*radius);
//			Y0 = BesselY0(WAVE*r);
//			Y0prime = BesselY0(WAVE*radius);
//			return ( J0prime*Y0 - Y0prime*J0 )/Besselnorm;
//			return 0.0;
			return BesselY0(WAVE*r)/BesselY0(WAVE*radius);
		}
		else
		{
//			return 0.0;
			return BesselY0(WAVE*r);

			J0 = BesselJ0(WAVE*r);
//	J0prime = BesselJ0(WAVE*radius);
			Y0 = BesselY0(WAVE*r);
//	Y0prime = BesselY0(WAVE*radius);
			return ( Y0*J0prime - J0*Y0prime )/Besselnorm;
	
//	return 0.0;
//	return BesselJ0(Lambda[0][1]*r);

//	return BesselJ0(Lambda[0][1]*r) + BesselJ0(Lambda2[0][1]*r2);
//	return BesselJ0(Lambda[0][1]*r) + BesselJ0(Lambda2[0][2]*r2) + BesselJ0(Lambda[0][3]*r);
//	return BesselJ0(Lambda[0][3]*r) + 2.0*BesselJ0(Lambda[0][3]*r2);
		}
	}
	else
	{
		if (INTERIOR_MODE==0)
		{
			for ( i = 1; i <= 3; i++)
				sum += ( pow((1.0/r),(double)(i)) * ( a[i-1]*cos(i*theta) + b[i-1]*sin(i*theta) ) );
			return sum;
		}
		else
		{
			for ( i = 1; i <= 4; i++)
				sum += ( pow(r,(double)(i)) * ( a[i-1]*cos(i*theta) + b[i-1]*sin(i*theta) ) );
			return sum;
		}
	}

////////////////////////////////////////////////////////////////////////

	for ( i = 1; i <= 4; i++)
		sum += pow(r,(double)(i)) * ( a[i-1]*cos(i*theta) + b[i-1]*sin(i*theta) );

	return sum;
}

double cal_partialu(double x, double y)
{
	double r, r2, temptheta, theta;
	double Dx, Dy, Dnorm, inner_pd, inner_pd2;

	if (KITE_TEST == 1)
	{
		cal_interpgrad(x,y,DIST_DX, DIST_DY, &Dx, &Dy);
		Dnorm = sqrt(Dx*Dx + Dy*Dy);
		Dx = Dx/Dnorm;
		Dy = Dy/Dnorm;

		inner_pd = x * PLANEWAVEX + y * PLANEWAVEY;
		inner_pd2 = -1. * (PLANEWAVEX * Dx + PLANEWAVEY * Dy);

//		return -1. * WAVE * inner_pd2 * sin(WAVE*inner_pd);
		return WAVE * inner_pd2 * sin(WAVE*inner_pd);
//		return -1. * WAVE * sin(WAVE*x) * Dx;
	}


	r = sqrt(pow(x-Cx,2.0) + pow(y-Cy,2.0));
	r2 = sqrt(pow(x-Cx2,2.0) + pow(y-Cy2,2.0));
	temptheta = atan((y-Cy)/(x-Cx));	// note the range of inverse trig functions is tricky	//
	if ((x >= Cx) && (y < Cy))	//	4th quadrant
		theta = temptheta + 2.0*PI;
	else if ((x >= Cx) && (y >= Cy))	//	1st quadrant
		theta = temptheta;
	else
		theta = temptheta + PI;

	if (LAPLACE == 1)
	{
		printf("Shouldn't need to use this program to solve Laplace Neumann Problem yet.\n");
		exit(0);
	}
	else
	{
		if (INTERIOR_MODE == 1)
		{
			return BesselJ1(WAVE*r)/BesselY0(WAVE*radius);	
		}
		else
		{
			return BesselJ1(WAVE*r) * WAVE;
		}
	}
	printf("Why here in cal_partialu?\n");
	exit(0);

}


double cal_partialuimg(double x, double y)
{
	double r, r2, temptheta, theta;
	double Dx, Dy, Dnorm, inner_pd, inner_pd2;

	if (KITE_TEST == 1)
	{
		cal_interpgrad(x,y,DIST_DX, DIST_DY, &Dx, &Dy);

		Dnorm = sqrt(Dx*Dx + Dy*Dy);
		Dx = Dx/Dnorm;
		Dy = Dy/Dnorm;

		inner_pd = x * PLANEWAVEX + y * PLANEWAVEY;
		inner_pd2 = -1. * (PLANEWAVEX * Dx + PLANEWAVEY * Dy);

//		return  WAVE * inner_pd2 * cos(WAVE*inner_pd);
		return  -1. * WAVE * inner_pd2 * cos(WAVE*inner_pd);
//		return -1. * WAVE * cos(WAVE*x) * Dx;
//		return WAVE * cos(WAVE*x) * Dx;
	}

	r = sqrt(pow(x-Cx,2.0) + pow(y-Cy,2.0));
	r2 = sqrt(pow(x-Cx2,2.0) + pow(y-Cy2,2.0));
	temptheta = atan((y-Cy)/(x-Cx));	// note the range of inverse trig functions is tricky	//
	if ((x >= Cx) && (y < Cy))	//	4th quadrant
		theta = temptheta + 2.0*PI;
	else if ((x >= Cx) && (y >= Cy))	//	1st quadrant
		theta = temptheta;
	else
		theta = temptheta + PI;

	if (LAPLACE == 1)
	{
		printf("Shouldn't need to use this program to solve Laplace Neumann Problem yet.\n");
		exit(0);
	}
	else
	{
		if (INTERIOR_MODE == 1)
		{
			return BesselY1(WAVE*r)/BesselY0(WAVE*radius);
		}
		else
		{
			return BesselY1(WAVE*r) * WAVE;
		}
	}
	printf("Why here in cal_partialuimg?\n");
	exit(0);
}


double cal_vdl(double x, double y)
{
	int i, j;
	static double a[4] = COEFFICIENTA;
	static double b[4] = COEFFICIENTB;
	static double Lambda[6][6] = ROOTBESSEL;
	static double Lambda2[6][6] = ROOTBESSEL;
	static int lambda_flag = 0;
	double r, r2, temptheta, theta;
	double sum = 0.0;
	double J0 = 0.0, J0prime = 0.0, Y0 = 0.0, Y0prime = 0.0;

//	return 1.0;

	printf("Deep in cal_vdl.\n");
	
	if (lambda_flag == 0)
	{
		for ( i = 0; i < 6; i++)
			for ( j = 0; j < 6; j++)
			{
				Lambda[i][j] /= radius;
				Lambda2[i][j] /= radius2;
			}
		lambda_flag++;
	}

	r = sqrt(pow(x-Cx,2.0) + pow(y-Cy,2.0));
	r2 = sqrt(pow(x-Cx2,2.0) + pow(y-Cy2,2.0));
	temptheta = atan((y-Cy)/(x-Cx));	// note the range of inverse trig functions is tricky	//
	if ((x >= Cx) && (y < Cy))	//	4th quadrant
		theta = temptheta + 2.0*PI;
	else if ((x >= Cx) && (y >= Cy))	//	1st quadrant
		theta = temptheta;
	else
		theta = temptheta + PI;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (LAPLACE == 0)
	{
		if (INTERIOR_MODE == 1)
		{
//			return 0.0;
			return BesselY0(WAVE*r)/BesselY0(WAVE*radius);
		}
		else
			return BesselY0(WAVE*r);
	}
	else
	{
		if (INTERIOR_MODE==0)
		{
			for ( i = 1; i <= 3; i++)
				sum += ( pow((r/(radius*radius)),(double)(i)) * ( a[i-1]*cos(i*theta) + b[i-1]*sin(i*theta) ) );
			return -1.0*sum;
		}
		else
		{
			for ( i = 1; i <= 4; i++)
				sum += ( pow(radius*radius/r,(double)(i)) * ( a[i-1]*cos(i*theta) + b[i-1]*sin(i*theta) ) );
			return -1.0*sum;	
		}
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for ( i = 1; i <= 4; i++)
		sum += pow(radius*radius/r,(double)(i)) * ( a[i-1]*cos(i*theta) + b[i-1]*sin(i*theta) );

	return -1.0*sum;
}

double cal_vdlimg(double x, double y)
{
	printf("cal_vdlimg not written yet.\n");
	return 0.0;
}


int NeumannTest(int size, double delta0, double epsilon0, double samplex, double sampley)
{
	int i, j;
	double zx, zy, zstarx, zstary, delta, thisJ, thisphi, thisphiimg, temp, thisnorm, result;
	double gradx, grady, gradnorm;
	double samplexstar, sampleystar, samplephi, samplephiimg, samplenorm, sampleu, sampleuimg, sampleunorm;
	double sampledist, samplegradx, samplegrady, samplegradnorm;
	double nearx, neary, neardist, nearu, nearuimg, nearphi, nearphiimg, nearnorm, nearunorm;
	double farx, fary, fardist, faru, faruimg, farphi, farphiimg, farnorm, farunorm;

	double fd_u, fd_uimg, fd_unorm;
	double tempvalue[2];
	double exactu, exactuimg, exactunorm, erroru, erroruimg, error, newerroru, newerroruimg, newerror;
	double exactnearu, exactnearuimg, exactnearunorm, exactfaru, exactfaruimg, exactfarunorm;
	double exactsampleu, exactsampleuimg, exactsampleunorm, sampleuerror, sampleuimgerror, sampleunormerror;

	int counti;
	double thisdist, thisgradx, thisgrady;
	static double tau, threshold;
//	epsilon = FACTOR * H;
	tau = H/TAU_FACTOR;
	threshold = tau * tau;

	exactsampleu = cal_u(samplex, sampley);
	exactsampleuimg = cal_uimg(samplex, sampley);
	exactsampleunorm = sqrt( exactsampleu*exactsampleu + exactsampleuimg*exactsampleuimg );

	if (BOUNDARY_UTEST == 1)
	{
		sampledist = cal_dist(samplex, sampley, DIST, DIST_DX, DIST_DY);
		cal_interpgrad(samplex, sampley, DIST_DX, DIST_DY, &samplegradx, &samplegrady);
		samplegradnorm = sqrt(samplegradx*samplegradx + samplegrady*samplegrady);
		samplegradx = samplegradx/samplegradnorm;
		samplegrady = samplegrady/samplegradnorm;
		samplexstar = samplex - sampledist * samplegradx;
		sampleystar = sampley - sampledist * samplegrady;


		exactu = cal_u(samplexstar, sampleystar);
		exactuimg = cal_uimg(samplexstar, sampleystar);
		exactunorm = sqrt( exactu*exactu + exactuimg*exactuimg );
	
		cal_interpgrad(samplexstar, sampleystar, DIST_DX, DIST_DY, &gradx, &grady);

		gradnorm = sqrt( gradx*gradx + grady*grady );

		printf("sample point (%lf, %lf), projected to (%lf, %lf)\n", samplex, sampley, samplexstar, sampleystar);


		if (ABS_FDTEST == 1)
		{
			neardist = ABS_NEARDIST;
			fardist = ABS_FARDIST;	
		}
		else
		{
			neardist = NEARFACTOR * H;
			fardist = FARFACTOR * H;
		}


		nearx = samplexstar - neardist * gradx/gradnorm;
		neary = sampleystar - neardist * grady/gradnorm;

		farx = samplexstar - fardist * gradx/gradnorm;
		fary = sampleystar - fardist * grady/gradnorm;

		//	printf("Near point (%lf, %lf), far point (%lf, %lf)\n", nearx, neary, farx, fary);
		exactnearu = cal_u( nearx, neary);
		exactnearuimg = cal_uimg( nearx, neary);
		exactnearunorm = sqrt( exactnearu*exactnearu + exactnearuimg*exactnearuimg );

		exactfaru = cal_u( farx, fary);
		exactfaruimg = cal_uimg( farx, fary);
		exactfarunorm = sqrt( exactfaru*exactfaru + exactfaruimg*exactfaruimg );
	}

	counti = 0;
	nearu = 0.;
	nearuimg = 0.;
	faru = 0.;
	faruimg = 0.;
	sampleu = 0.;
	sampleuimg = 0.;

//	if (ABS_MODE != 1)
//		epsilon = FACTOR * h;

	for ( i = 0; i < SIZE; i++)
	{
		zy = INT_L + H * i;
		for ( j = 0; j < SIZE; j++)		
		{
			zx = INT_L + H*j;

//			if (ABS_MODE == 1)
//			{
//				epsilon = ABS_EPSILON;
//			}
//			else if (ABS_MODE == 0)
//			{
//				thisgradx = gradientdx(zx, zy);
//				thisgrady = gradientdy(zx, zy);
//				epsilon = ( fabs(thisgradx) + fabs(thisgrady) ) * FACTOR * h;
//			}

//			thisdist = distance(zx, zy);
			thisdist = DIST[i][j];
//			printf("i = %d, j = %d, point (%lf, %lf), distance %lf.\n", i, j, zx, zy, thisdist);

			if (TRAD_ONESIDED == 0)
			{
				if (fabs(thisdist) > epsilon0)
					continue;
			}
			else
			{
				if ( (thisdist > epsilon0) || (thisdist < delta0) )
					continue;
			}

//			zstarx = starx(zx, zy);
//			zstary = stary(zx, zy);

			zstarx = zx - thisdist * DIST_DX[i][j];
			zstary = zy - thisdist * DIST_DY[i][j];

			if (TRAD_ONESIDED == 0)
				delta = cal_delta(epsilon0, thisdist);
			else
				delta = 2. * cal_delta(epsilon0, thisdist);

//			thisJ = cal_J(J_M, zx, zy);
			thisJ = 1. + thisdist * CURV[i][j];
			result = HSQUARED * delta * thisJ;

			if (BOUNDARY_UTEST == 1)
			{
				nearnorm = sqrt( (zstarx-nearx)*(zstarx-nearx) + (zstary-neary)*(zstary-neary) );
				if (nearnorm > tau)
				{
					nearphi = phi( nearx, neary, zstarx, zstary);
					nearphiimg = phiimg( nearx, neary, zstarx, zstary);
				}
				else
				{
					nearphi = 0.;
					nearphiimg = 0.;
				}
				nearu += ( ( nearphi * Density[counti] - nearphiimg * Densityimg[counti] ) * result );
				nearuimg += ( ( nearphi * Densityimg[counti] + nearphiimg * Density[counti] ) * result );


				farnorm = sqrt( (zstarx-farx)*(zstarx-farx) + (zstary-fary)*(zstary-fary) );
				if (farnorm > tau)
				{
					farphi = phi( farx, fary, zstarx, zstary);
					farphiimg = phiimg( farx, fary, zstarx, zstary);
				}
				else
				{
					farphi = 0.;
					farphiimg = 0.;
				}
				faru += ( ( farphi * Density[counti] - farphiimg * Densityimg[counti] ) * result );
				faruimg += ( ( farphi * Densityimg[counti] + farphiimg * Density[counti] ) * result );
			}
			samplenorm = sqrt( (zstarx-samplex)*(zstarx-samplex) + (zstary-sampley)*(zstary-sampley) );
			if (samplenorm > tau)
			{
				samplephi = BesselY0(WAVE * samplenorm)/(4.);
				samplephiimg = BesselJ0(WAVE * samplenorm)/(-4.);
			}
			else
			{
				samplephi = 0.;
				samplephiimg = 0.;
			}

			sampleu += ( (samplephi * Density[counti] - samplephiimg * Densityimg[counti]) * result );
			sampleuimg += ( (samplephi * Densityimg[counti] + samplephiimg * Density[counti]) * result );

			if ( (nearu != nearu) || (faru != faru) || (sampleu != sampleu) )
			{
				printf("At %d, Density = (%lf, %lf), nearphi (%lf, %lf), farphi(%lf, %lf), result %lf\n", \
				counti, Density[counti], Densityimg[counti], nearphi, nearphiimg, farphi, farphiimg, \
				result);
				printf("sample point (%lf, %lf), phi (%lf, %lf), u = %lf, uimg = %lf.\n", \
				samplex, sampley, samplephi, samplephiimg, sampleu, sampleuimg);
				exit(0);
			}

			counti++;
		}
	}

	if (BOUNDARY_UTEST == 1)
	{
		nearunorm = sqrt(nearu*nearu + nearuimg*nearuimg);
		farunorm = sqrt(faru*faru + faruimg*faruimg);
		fd_u = nearu + neardist * (nearu - faru)/(fardist-neardist);
		fd_uimg = nearuimg + neardist * (nearuimg - faruimg)/(fardist-neardist);

		fd_unorm = sqrt(fd_u*fd_u + fd_uimg*fd_uimg);

		erroru = fabs(fd_u - exactu)/fabs(exactu);
		erroruimg = fabs(fd_uimg - exactuimg)/fabs(exactuimg);
		error = fabs(fd_unorm - exactunorm)/fabs(exactunorm);
//	printf("FarR: exact = %lf, appr = %lf, error %.2E\n", exactfaru, faru, fabs(exactfaru-faru)/fabs(exactfaru));
//	printf("FaIm: exact = %lf, appr = %lf, error %.2E\n", exactfaruimg, faruimg, \
//								fabs(exactfaruimg-faruimg)/fabs(exactfaruimg));
//	printf("FarN: exact = %lf, appr = %lf, error %.2E\n", exactfarunorm, farunorm, \
//							fabs(exactfarunorm-farunorm)/fabs(exactfarunorm));
//
//	printf("NeRe: exact = %lf, appr = %lf, error %.2E\n", exactnearu, nearu, fabs(exactnearu-nearu)/fabs(exactnearu));
//	printf("NeIm: exact = %lf, appr = %lf, error %.2E\n", exactnearuimg, nearuimg, \
//								fabs(exactnearuimg-nearuimg)/fabs(exactnearuimg));
//	printf("NeaN: exact = %lf, appr = %lf, error %.2E\n", exactnearunorm, nearunorm, \
//					fabs(exactnearunorm-nearunorm)/fabs(exactnearunorm));

//	printf("Real: exact = %lf, fd_u = %lf, error %.2E\n", exactu, fd_u, erroru);
//	printf("Imag: exact = %lf, fd_u = %lf, error %.2E\n", exactuimg, fd_uimg, erroruimg);
//	printf("Norm: exact = %lf, fd_u = %lf, error %.2E\n", exactunorm, fd_unorm, error);
	}

	sampleunorm = sqrt(sampleu*sampleu + sampleuimg*sampleuimg);
	sampleuerror = fabs(sampleu - exactsampleu)/fabs(exactsampleu);
	sampleuimgerror = fabs(sampleuimg - exactsampleuimg)/fabs(exactsampleuimg);
	sampleunormerror = fabs(sampleunorm - exactsampleunorm)/fabs(exactsampleunorm);
	
	printf("Real: exact = %lf, trad_approxu = %lf, error %.2E\n", exactsampleu, sampleu, sampleuerror);
	printf("Imag: exact = %lf, trad_approxu = %lf, error %.2E\n", exactsampleuimg, sampleuimg, sampleuimgerror);
	printf("Norm: exact = %lf, trad_approxu = %lf, error %.2E\n", exactsampleunorm, sampleunorm, sampleunormerror);
	return 0;
}


int New_NeumannTest(int size, double samplex, double sampley)
{
	int i, j;
	double zx, zy, zstarx, zstary, delta, thisJ, thisphi, thisphiimg, temp, thisnorm, result;
	double samplexstar, sampleystar, samplegradx, samplegrady, sampledist, samplegradnorm;
	double *newx1, *newx2, *xminusy1, *xminusy2, *samplekernel, *samplekernelimg;
	double threshold = 0.01*H*H;
	double polysum, polysumimg, sinesum, sinesumimg, polyu, polyuimg, polyunorm, sineu, sineuimg, sineunorm;
	double polyerroru, polyerroruimg, polyerrorunorm, sineerroru, sineerroruimg, sineerrorunorm;
	double exactu, exactuimg, exactunorm;

	double tau;

	tau = H/TAU_FACTOR;
	threshold = tau * tau;


	newx1 = (double *) malloc(size * sizeof(double));
	newx2 = (double *) malloc(size * sizeof(double));
	xminusy1 = (double *) malloc(size * sizeof(double));
	xminusy2 = (double *) malloc(size * sizeof(double));
	samplekernel = (double *) malloc(size * sizeof(double));
	samplekernelimg = (double *) malloc(size * sizeof(double));

	sampledist = cal_dist(samplex, sampley, DIST, DIST_DX, DIST_DY);
	cal_interpgrad(samplex, sampley, DIST_DX, DIST_DY, &samplegradx, &samplegrady);
	samplegradnorm = sqrt(samplegradx*samplegradx + samplegrady*samplegrady);
	samplegradx = samplegradx/samplegradnorm;
	samplegrady = samplegrady/samplegradnorm;
	samplexstar = samplex - sampledist * samplegradx;
	sampleystar = sampley - sampledist * samplegrady;
//	samplexstar = starx(samplex, sampley);
//	sampleystar = stary(samplex, sampley);

//	samplegradx = cal_interpgrad(1, samplexstar, sampleystar);
//	samplegrady = cal_interpgrad(2, samplexstar, sampleystar);
	cal_interpgrad(samplexstar, sampleystar, DIST_DX, DIST_DY, &samplegradx, &samplegrady);

	exactu = cal_u(samplexstar, sampleystar);
	exactuimg = cal_uimg(samplexstar, sampleystar);
	exactunorm = sqrt( exactu*exactu + exactuimg*exactuimg );

	for ( i = 0; i < size; i++)
	{
		newx1[i] = samplexstar - samplegradx * fabs(newdist[i]);
		newx2[i] = sampleystar - samplegrady * fabs(newdist[i]);

		xminusy1[i] = newx1[i] - newzstarx[i];
		xminusy2[i] = newx2[i] - newzstary[i];

		thisnorm = sqrt( xminusy1[i]*xminusy1[i] + xminusy2[i]*xminusy2[i] );
		if (thisnorm > threshold)
		{
			samplekernel[i] = BesselY0(WAVE * thisnorm)/4.;
			samplekernelimg[i] = BesselJ0(WAVE * thisnorm)/(-4.);
		}
		else
		{
			samplekernel[i] = 0.;
			samplekernelimg[i] = 0.;
		}
	}

	polysum = 0.;
	polysumimg = 0.;
	sinesum = 0.;
	sinesumimg = 0.;
	for ( i = 0; i < size; i++)
	{
		polysum +=  (( samplekernel[i] * newPolyDensity[i] - samplekernelimg[i] * newPolyDensityimg[i] )  * \
			PolyResult[i] );
		polysumimg += ( ( samplekernelimg[i] * newPolyDensity[i] + samplekernel[i] * newPolyDensityimg[i] ) * \
			PolyResult[i] );
		sinesum += ( ( samplekernel[i] * newSineDensity[i] - samplekernelimg[i] * newSineDensityimg[i] ) * \
			SineResult[i] );
		sinesumimg += ( ( samplekernelimg[i] * newSineDensity[i] + samplekernel[i] * newSineDensityimg[i] ) * \
			SineResult[i] );
//		printf("%d, Den %lf, Denimg %lf, PolyResult %lf\n",\
//			 i, newPolyDensity[i], newPolyDensityimg[i], PolyResult[i]);
	}
	polyu = polysum;
	polyuimg = polysumimg;
	polyunorm = sqrt( polyu * polyu + polyuimg * polyuimg );
	sineu = sinesum;
	sineuimg = sinesumimg;
	sineunorm = sqrt( sineu * sineu + sineuimg * sineuimg );

	polyerroru = fabs(polyu - exactu)/fabs(exactu);
	polyerroruimg = fabs(polyuimg - exactuimg)/fabs(exactuimg);
	polyerrorunorm = fabs(polyunorm - exactunorm)/fabs(exactunorm);

	sineerroru = fabs(sineu - exactu)/fabs(exactu);
	sineerroruimg = fabs(sineuimg - exactuimg)/fabs(exactuimg);
	sineerrorunorm = fabs(sineunorm - exactunorm)/fabs(exactunorm);
	
//	printf("Poly:\n");
	printf("Real: exact = %lf, polyu %lf, error %.2E, sineu %lf, error %.2E\n", \
	exactu, polyu, polyerroru, sineu, sineerroru);
	printf("Imag: exact = %lf, polyu %lf, error %.2E, sineu %lf, error %.2E\n", \
	exactuimg, polyuimg, polyerroruimg, sineuimg, sineerroruimg);
	printf("Norm: exact = %lf, polyu %lf, error %.2E, sineu %lf, error %.2E\n", \
	exactunorm, polyunorm, polyerrorunorm, sineunorm, sineerrorunorm);

	free(newx1);
	free(newx2);
	free(xminusy1);
	free(xminusy2);
	free(samplekernel);
	free(samplekernelimg);

	return 0;

}


///	EVALUATE U AT A POINT OFF OR ON INTERFACE.
//	IF A POINT IS ON INTERFACE, BOTH DOUBLE AND SINGLE LAYER USE THE DYNAMIC ESTIMATOR (X CHANGES WITH Y IN INTEGRAL)
//	IF A POINT IS OFF INTERFACE, DO NOT USE DYNAMIC ESTIMATOR (X STAYS FIXED FOR ANY Y IN INTEGRAL)

int Combo_NeumannTest(int size, double samplex, double sampley, double eta)
{
	int i, j;

	double samplexstar, sampleystar, samplegradx, samplegrady, samplegradnorm, sampledist;
	double *newx1, *newx2, *xminusy1, *xminusy2, *xminusynorm, *samplekernel, *samplekernelimg;

	double *farfieldkernel, *farfieldkernelimg;
	double directionx1, directionx2, pointxnorm, constterm, polyfarfield, polyfarfieldimg, sinefarfield, sinefarfieldimg;

	double threshold = 0.01*H*H, tau;
	double polysum, polysumimg, sinesum, sinesumimg, polyu, polyuimg, polyunorm, sineu, sineuimg, sineunorm;
	double polyerroru, polyerroruimg, polyerrorunorm, sineerroru, sineerroruimg, sineerrorunorm;
	double exactu, exactuimg, exactunorm, exactustar, exactustarimg, exactustarnorm;
	double diffu, diffuimg, diffunorm, sanitysineu, sanitysineuimg, sanitypolyu, sanitypolyuimg;



	tau = H/TAU_FACTOR;
	threshold = tau* tau;	

	newx1 = (double *) malloc(size * sizeof(double));
	newx2 = (double *) malloc(size * sizeof(double));
	xminusy1 = (double *) malloc(size * sizeof(double));
	xminusy2 = (double *) malloc(size * sizeof(double));
	samplekernel = (double *) malloc(size * sizeof(double));
	samplekernelimg = (double *) malloc(size * sizeof(double));

	if (KITE_TEST == 1)
	{
		pointxnorm = sqrt(samplex*samplex + sampley*sampley);
		directionx1 = samplex/pointxnorm;
		directionx2 = sampley/pointxnorm;
		constterm = 0.25/ sqrt(WAVE*PI);

		farfieldkernel = (double *) malloc(size * sizeof(double));
		farfieldkernelimg = (double *) malloc(size * sizeof(double));
	}

	if (BOUNDARY_UTEST == 1)
	{
		sampledist = cal_dist(samplex, sampley, DIST, DIST_DX, DIST_DY);
		cal_interpgrad(samplex, sampley, DIST_DX, DIST_DY, &samplegradx, &samplegrady);
		samplegradnorm = sqrt(samplegradx*samplegradx + samplegrady*samplegrady);
		samplegradx = samplegradx/samplegradnorm;
		samplegrady = samplegrady/samplegradnorm;
		samplexstar = samplex - sampledist * samplegradx;
		sampleystar = sampley - sampledist * samplegrady;
		exactustar = cal_u(samplexstar, sampleystar);
		exactustarimg = cal_uimg(samplexstar, sampleystar);
		exactustarnorm = sqrt( exactustar*exactustar + exactustarimg*exactustarimg );
	}
	
	exactu = cal_u(samplex, sampley);
	exactuimg = cal_uimg(samplex, sampley);
	exactunorm = sqrt( exactu*exactu + exactuimg*exactuimg );

	//	OFF INTERFACE, u(x) = [ dG/dny(x,y*) + i * eta * G(x,y*) ] W(y) J(y) beta(y) H^2	//
	#pragma omp parallel for private(i) schedule(static, 16)
	for ( i = 0; i < size; i++)
	{
		if (KITE_TEST == 1)
		{
			double xdotny, xdoty, pointynorm, directiony1, directiony2;
			double kiteparameter, kiteSine, kiteCosine;

			xdotny = -1. * (directionx1 * newgradx[i] + directionx2 * newgrady[i]);
//			xdotny = (directionx1 * newgradx[i] + directionx2 * newgrady[i]);
//			pointynorm = sqrt( newzstarx[i] * newzstarx[i] + newzstary[i] * newzstary[i] );
//			directiony1 = newzstarx[i] / pointynorm;
//			directiony2 = newzstary[i] / pointynorm;

			xdoty = directionx1 * newzstarx[i] + directionx2 * newzstary[i];
			kiteparameter = WAVE * xdoty;
			kiteSine = sin(kiteparameter);
			kiteCosine = cos(kiteparameter);

			farfieldkernel[i] = constterm * (WAVE * xdotny - eta) * ( kiteCosine - kiteSine );
			farfieldkernelimg[i] = -1. * constterm * (WAVE * xdotny - eta) * ( kiteSine + kiteCosine );
//			farfieldkernel[i] = constterm * (WAVE * xdotny + eta) * ( kiteCosine - kiteSine );
//			farfieldkernelimg[i] = -1. * constterm * (WAVE * xdotny + eta) * ( kiteSine + kiteCosine );
//			if (i < 10)
//			{
//				printf("x (%lf, %lf) grad (%lf, %lf), y (%lf, %lf) grad(%lf, %lf)\n", \
//				samplex, sampley, directionx1, directionx2, newzstarx[i], newzstary[i], directiony1, \
//				directiony2);
//
//				printf("constterm %lf, xdoty %lf, xdotny %lf, para %lf, kitesine %lf, kitecosine %lf\n", \
//				constterm, xdoty, xdotny, kiteparameter, kiteSine, kiteCosine);
//			}
		}
		double xmy1, xmy2, thisnorm;
		double partial, partialimg;
		double commonterm1, commonterm2, parameter;
		double singlereal, singleimg;

		xmy1 = samplex - newzstarx[i];
		xmy2 = sampley - newzstary[i]; 
		

		thisnorm = sqrt( xmy1*xmy1 + xmy2*xmy2 );
		parameter = WAVE * thisnorm;

		//	k Y1(k|x-y*|) / 4|x-y*|	//
		commonterm1 = WAVE * BesselY1(parameter) / (4.*thisnorm);
		//	k J1(k|x-y*|) / 4|x-y*|	//
		commonterm2 = WAVE * BesselJ1(parameter) / (-4.*thisnorm);

		//	dG/dny = -1 * Grad G dot grad d = -ik/4 * H1(k|x-y*|)/|x-y*| * [(x-y) dot grad d]	//
		//	There's a negative from H0' = -H1, one from ny = -grad d, one from d|x-y*|/dy* = -(x-y*)/|x-y*|
		partial = -1. * commonterm1 * ( xmy1 * newgradx[i] + xmy2 * newgrady[i] );
		partialimg = -1. * commonterm2 * ( xmy1 * newgradx[i] + xmy2 * newgrady[i] );


		//	-i * eta * G(x,y) = -i eta * (-i/4)*(J0 + iY0) = 0.25*eta*(J0+iY0)	//
		singlereal = 0.25 * eta * BesselJ0(parameter);
		singleimg = 0.25 * eta * BesselY0(parameter);
		
		samplekernel[i] = -1.*(partial + singlereal);
		samplekernelimg[i] = -1.*(partialimg + singleimg);

	}

	polysum = 0.0;
	polysumimg = 0.0;
	sinesum = 0.0;
	sinesumimg = 0.0;
	polyfarfield = 0.0;
	polyfarfieldimg = 0.0;
	sinefarfield = 0.0;
	sinefarfieldimg = 0.0;

	#pragma omp parallel for private(i) schedule(static, 16) \
		reduction(+: polysum, polysumimg, polyfarfield, polyfarfieldimg, \
				sinesum, sinesumimg, sinefarfield, sinefarfieldimg)
	for ( i = 0 ; i < size; i++)
	{
		if (POLY_TEST == 1)
		{
		polysum += ( (samplekernel[i] * newPolyDensity[i] - samplekernelimg[i] * newPolyDensityimg[i] ) * \
				PolyResult[i] );
		polysumimg += ( (samplekernelimg[i] * newPolyDensity[i] + samplekernel[i] * newPolyDensityimg[i]) * \
				PolyResult[i] );
		if (KITE_TEST == 1)
		{
			polyfarfield += ( (farfieldkernel[i]*newPolyDensity[i]-farfieldkernelimg[i]*newPolyDensityimg[i]) *\
					PolyResult[i] );
			polyfarfieldimg += ( (farfieldkernelimg[i]*newPolyDensity[i] + \
					farfieldkernel[i]*newPolyDensityimg[i]) * PolyResult[i] );
		}
		}
		if (SINE_TEST == 1)
		{
		sinesum += ( (samplekernel[i] * newSineDensity[i] - samplekernelimg[i] * newSineDensityimg[i] ) * \
				SineResult[i] );
		sinesumimg += ( (samplekernelimg[i] * newSineDensity[i] + samplekernel[i] * newSineDensityimg[i]) * \
				SineResult[i] );
		if (KITE_TEST == 1)
		{
			sinefarfield += ( (farfieldkernel[i]*newSineDensity[i]-farfieldkernelimg[i]*newSineDensityimg[i]) *\
					SineResult[i] );
			sinefarfieldimg += ( (farfieldkernelimg[i]*newSineDensity[i] + \
					farfieldkernel[i]*newSineDensityimg[i]) * SineResult[i] );
		}
		}
	}

	if (POLY_TEST == 1)
	{
		polyu = polysum;
		polyuimg = polysumimg;
		polyunorm = sqrt( polyu * polyu + polyuimg * polyuimg );
		polyerroru = fabs(polyu - exactu)/fabs(exactu);
		polyerroruimg = fabs(polyuimg - exactuimg)/fabs(exactuimg);
		polyerrorunorm = fabs(polyunorm - exactunorm)/fabs(exactunorm);
	}
	if (SINE_TEST == 1)
	{
		sineu = sinesum;
		sineuimg = sinesumimg;
		sineunorm = sqrt( sineu * sineu + sineuimg * sineuimg );
		sineerroru = fabs(sineu - exactu)/fabs(exactu);
		sineerroruimg = fabs(sineuimg - exactuimg)/fabs(exactuimg);
		sineerrorunorm = fabs(sineunorm - exactunorm)/fabs(exactunorm);
	}
	
//	printf("Poly:\n");
	printf("Real: exact = %lf, ", exactu);
	if (POLY_TEST == 1)
		printf("polyu %lf, error %.2E, ", polyu, polyerroru);
	if (SINE_TEST == 1)
		printf("sineu %lf, error %.2E\n", sineu, sineerroru);
	printf("Imag: exact = %lf, ", exactuimg);
	if (POLY_TEST == 1)
		printf("polyu %lf, error %.2E, ", polyuimg, polyerroruimg);
	if (SINE_TEST == 1)
		printf("sineu %lf, error %.2E\n", sineuimg, sineerroruimg);
	printf("Norm: exact = %lf, ", exactunorm);
	if (POLY_TEST == 1)
		printf("polyu %lf, error %.2E, ", polyunorm, polyerrorunorm);
	if (SINE_TEST == 1)
		printf("sineu %lf, error %.2E\n", sineunorm, sineerrorunorm);

	if (KITE_TEST == 1)
	{
		printf("Farfield u at direction (%lf, %lf): ", directionx1, directionx2);
		if (POLY_TEST == 1)
		{
			printf("poly %lf + %lfi, ", polyfarfield, polyfarfieldimg);
		}
		if (SINE_TEST == 1)
		{
			printf("sine %lf + %lfi\n", sinefarfield, sinefarfieldimg);

			sanitysineu = ( cos(WAVE*pointxnorm) * sinefarfield - sin(WAVE*pointxnorm) * sinefarfieldimg)/ \
					sqrt(pointxnorm);
			sanitysineuimg = ( cos(WAVE*pointxnorm) * sinefarfieldimg + sin(WAVE*pointxnorm) * sinefarfield)/ \
					sqrt(pointxnorm);

			diffu = sineu - sanitysineu;
			diffuimg = sineuimg - sanitysineuimg;
			diffunorm = sqrt(diffu*diffu + diffuimg*diffuimg);

			printf("Sanity check u difference %lf + %lfi, norm %lf.\n", diffu, diffuimg, diffunorm);
					
		}

		free(farfieldkernel);
		free(farfieldkernelimg);
	}




	free(newx1);
	free(newx2);
	free(xminusy1);
	free(xminusy2);
	free(samplekernel);
	free(samplekernelimg);

	return 0;
		

/*
		
		//	nx dot ny	//
		inner_pd1 = samplegradx * newgradx + samplegrady * newgrady;
		//	nx dot (x-y*)	//
		inner_pd2 = samplegradx * xmy1 + samplegrady * xmy2;
		//	ny dot (x-y*)	//
		inner_pd3 = xmy1 * newgradx + xmy2 * newgrady;

		//	k|x-y*|		//
		parameter = WAVE * thisnorm;
		commonterm1 = inner_pd1 / thisnorm;
		commonterm2 = inner_pd2 * inner_pd3 / (thisnorm * thisnorm);
		BesselY1value = BesselY1( parameter );
		BesselJ1value = BesselJ1( parameter );
		//	first term
	}

	for (int i = 0; i < size; i++)
	{
		newx1[i] = samplexstar - samplegradx * fabs(newdist[i]);
		newx2[i] = sampleystar - samplegrady * fabs(newdist[i]);

		xminusy1[i] = newx1[i] - newzstarx[i];
		xminusy2[i] = newx2[i] - newzstary[i];
		thisnorm = sqrt( xminusy1[i]*xminusy1[i] + xminusy2[i]*xminusy2[i] );

		//	nx dot ny	//
		inner_pd1 = samplegradx * newgradx[i] + samplegrady * newgrady[i];
		//	(x-y*) dot nx	//
		inner_pd2 = samplegradx * xminusy1[i] + samplegrady * xminusy2[i];
		//	(x-y*) dot ny	//
		inner_pd3 = newgradx[i] * xminusy1[i] + newgrady[i] * xminusy2[i];

		parameter = thisnorm * WAVE;
		commonterm1 = inner_pd1 / thisnorm;
		commonterm2 = inner_pd2 * inner_pd3 / (thisnorm * thisnorm);
		BesselY1value = BesselY1( parameter );
		BesselJ1value = BesselJ1( parameter );

		//	first term, H1(k|x-y*|) * [nx dot ny] / |x-y*|	(divided by k/2) 	//
		term1R = BesselY1value * commonterm1;
		term1I = BesselJ1value * commonterm1;
		//	second term, k * H0(k|x-y*|) * [(x-y*) dot nx] * [(x-y*) dot ny]/ |x-y*|^2	//
		term2R = BesselY0(parameter) * WAVE * commonterm2;
		term2I = BesselJ0(parameter) * WAVE * commonterm2;
		//	third term, H1(k|x-y*|) * [(x-y*) dot nx] * [(x-y*) dot ny] / |x-y*|^3
		term3R = BesselY1value * commonterm2 / xminusynorm;
		term3I = BesselJ1value * commonterm2 / xminusynorm;

		


		if (thisnorm > threshold)
		{
			singlereal = BesselY0(WAVE * thisnorm)/4.;
			singleimg = BesselY0(WAVE * thisnorm)/(-4.);

			doublereal = WAVE * (term1R + term2R + term3R)/2.;
			doubleimg = WAVE * (term1I + term2I + term3I)/(-2.);
			

			samplekernel[i] = BesselY0(WAVE * thisnorm)/4.;
			samplekernelimg[i] = BesselJ0(WAVE * thisnorm)/(-4.);
		}
		else
		{
			samplekernel[i] = 0.;
			samplekernelimg[i] = 0.;
		}
	}

	polysum = 0.;
	polysumimg = 0.;
	sinesum = 0.;
	sinesumimg = 0.;
	for (int i = 0; i < size; i++)
	{
		polysum +=  (( samplekernel[i] * newPolyDensity[i] - samplekernelimg[i] * newPolyDensityimg[i] )  * \
			PolyResult[i] );
		polysumimg += ( ( samplekernelimg[i] * newPolyDensity[i] + samplekernel[i] * newPolyDensityimg[i] ) * \
			PolyResult[i] );
		sinesum += ( ( samplekernel[i] * newSineDensity[i] - samplekernelimg[i] * newSineDensityimg[i] ) * \
			SineResult[i] );
		sinesumimg += ( ( samplekernelimg[i] * newSineDensity[i] + samplekernel[i] * newSineDensityimg[i] ) * \
			SineResult[i] );
//		printf("%d, Den %lf, Denimg %lf, PolyResult %lf\n",\
//			 i, newPolyDensity[i], newPolyDensityimg[i], PolyResult[i]);
	}
	polyu = polysum;
	polyuimg = polysumimg;
	polyunorm = sqrt( polyu * polyu + polyuimg * polyuimg );
	sineu = sinesum;
	sineuimg = sinesumimg;
	sineunorm = sqrt( sineu * sineu + sineuimg * sineuimg );

	polyerroru = fabs(polyu - exactu)/fabs(exactu);
	polyerroruimg = fabs(polyuimg - exactuimg)/fabs(exactuimg);
	polyerrorunorm = fabs(polyunorm - exactunorm)/fabs(exactunorm);

	sineerroru = fabs(sineu - exactu)/fabs(exactu);
	sineerroruimg = fabs(sineuimg - exactuimg)/fabs(exactuimg);
	sineerrorunorm = fabs(sineunorm - exactunorm)/fabs(exactunorm);
	
//	printf("Poly:\n");
	printf("Real: exact = %lf, polyu %lf, error %.2E, sineu %lf, error %.2E\n", \
	exactu, polyu, polyerroru, sineu, sineerroru);
	printf("Imag: exact = %lf, polyu %lf, error %.2E, sineu %lf, error %.2E\n", \
	exactuimg, polyuimg, polyerroruimg, sineuimg, sineerroruimg);
	printf("Norm: exact = %lf, polyu %lf, error %.2E, sineu %lf, error %.2E\n", \
	exactunorm, polyunorm, polyerrorunorm, sineunorm, sineerrorunorm);

	free(newx1);
	free(newx2);
	free(xminusy1);
	free(xminusy2);
	free(samplekernel);
	free(samplekernelimg);

	return 0;
*/

}


int print_SINGLE_kernel(int size, double samplex, double sampley)
{
	int i, j;
	FILE *fpk, *fpkreg, *fpkimg, *fpkimgreg;
	int counti = 0;
	double zx, zy, zstarx, zstary, delta;
//	double partial, temp, partialimg;
//	double partialx, partialy, gradx, grady, partialximg, partialyimg;
//	double partial_reg, partialimg_reg;
	double thisJ, thisphi, thisphiimg, temp, thisnorm;
	double tempvalue[2];
	double thisdist, thisgradx, thisgrady;

	static double epsilon, tau, threshold;
	epsilon = FACTOR * H;
	tau = H/TAU_FACTOR;
	threshold = tau * tau;

	if ((fpk = fopen("largesamplekernel.txt", "w+"))==NULL)
	{
		printf("Open file sample.txt error.\n");
		return 0;
	}
	if ((fpkreg = fopen("largesamplekernelreg.txt", "w+"))==NULL)
	{
		printf("Open file samplekernelreg.txt error.\n");
		return 0;
	}
	if ((fpkimg = fopen("largesamplekernelimg.txt", "w+"))==NULL)
	{
		printf("Open file samplekernelimg.txt error.\n");
		return 0;
	}
	if ((fpkimgreg = fopen("largesamplekernelimgreg.txt", "w+"))==NULL)
	{
		printf("Open file samplekernelimgreg.txt error.\n");
		return 0;
	}
	fprintf(fpkimg, "%d\n", size);
	fprintf(fpkimgreg, "%d\n", size);
	fprintf(fpkimg, "%.15lf %.15lf\n", samplex, sampley);
	fprintf(fpkimgreg, "%.15lf %.15lf\n", samplex, sampley);

	if (cal_dist(samplex, sampley, DIST, DIST_DX, DIST_DY) <= 0.0)	//	OUTSIDE		//
	{
		fprintf(fpkimg, "%.15lf\n", cal_uimg(samplex, sampley));
		fprintf(fpkimgreg, "%.15lf\n", cal_uimg(samplex, sampley));
	}
	else					//	INSIDE		//
	{
		fprintf(fpkimg, "%.15lf\n", cal_vdlimg(samplex, sampley));
		fprintf(fpkimgreg, "%.15lf\n", cal_vdlimg(samplex, sampley));
	}

	fprintf(fpk, "%d\n", size);
	fprintf(fpkreg, "%d\n", size);
	fprintf(fpk, "%.15lf %.15lf\n", samplex, sampley);
	fprintf(fpkreg, "%.15lf %.15lf\n", samplex, sampley);

	if (cal_dist(samplex, sampley, DIST, DIST_DX, DIST_DY) <= 0.0)	//	OUTSIDE		//
	{
		fprintf(fpk, "%.15lf\n", cal_u(samplex, sampley));
		fprintf(fpkreg, "%.15lf\n", cal_u(samplex, sampley));
	}
	else					//	INSIDE		//
	{
		fprintf(fpk, "%.15lf\n", cal_vdl(samplex, sampley));
		fprintf(fpkreg, "%.15lf\n", cal_vdl(samplex, sampley));
	}


	if (ABS_MODE != 1)
		epsilon = FACTOR * H;

	for (i = 0; i < SIZE; i++)
	{
		zy = INT_L + H * i;
		for ( j = 0; j < SIZE; j++)		
		{
			zx = INT_L + H*j;
				
			if (ABS_MODE == 1)
			{
				epsilon = ABS_EPSILON;
			}
			else if (ABS_MODE == 0)
			{
				thisgradx = gradientdx(zx, zy);
				thisgrady = gradientdy(zx, zy);
//				epsilon = ( fabs(thisgradx) + fabs(thisgrady) ) * FACTOR * H;
			}

//			thisdist = distance(zx, zy);
			thisdist = DIST[i][j];
//			if (fabs(thisdist) > epsilon)
			if (( thisdist > epsilon) || (thisdist < 0.) )
				continue;
			counti++;
//			zstarx = starx(zx, zy);
//			zstary = stary(zx, zy);

			zstarx = zx - thisdist - DIST_DX[i][j];
			zstary = zy - thisdist - DIST_DY[i][j];

			delta = cal_delta(epsilon, thisdist);
//			thisJ = cal_J(J_M, zx, zy);
			thisJ = 1. + thisdist * CURV[i][j];

			thisnorm =  sqrt( (zstarx-samplex)*(zstarx-samplex) + (zstary-sampley)*(zstary-sampley) );

			if (thisnorm > tau)
			{
				thisphi = phi(samplex, sampley, zstarx, zstary);
				thisphiimg = phiimg(samplex, sampley, zstarx, zstary);
			}
			else
			{
				thisphi = 0.;
				thisphiimg = 0.;
			}

		

			temp = H*H * thisphi * delta * thisJ;
			fprintf(fpk, "%.15lf ", temp);
			temp = H*H * thisphi * delta * thisJ;
			fprintf(fpkreg, "%.15lf ", temp);
//			if (LAPLACE == 0)
//			{
			temp = H*H * thisphiimg * delta * thisJ;
			fprintf(fpkimg, "%.15lf ", temp);
			temp = H*H * thisphiimg * delta * thisJ;
			fprintf(fpkimgreg, "%.15lf ", temp);
//			}
		}
	}
	fclose(fpk);
	fclose(fpkreg);
	fclose(fpkimg);
	fclose(fpkimgreg);

	return 0;
}


int print_kernel(int size, double samplex, double sampley)
{
	int i, j;
	FILE *fpk, *fpkreg, *fpkimg, *fpkimgreg;
	int counti = 0;
	double zx, zy, zstarx, zstary, delta, thisJ;
	double partial, temp, partialimg;
	double partialx, partialy, gradx, grady, partialximg, partialyimg;
	double partial_reg, partialimg_reg;
	double tempvalue[2];
	double thisdist, thisgradx, thisgrady;

	static double epsilon, tau, threshold, regfactor, regfactorimg;

	epsilon = FACTOR * H;
	tau = H/TAU_FACTOR;
	threshold = tau * tau;
	regfactor = 1.0/(4.0*PI*radius);
	regfactorimg = 0.0;

	

	if ((fpk = fopen("largesamplekernel.txt", "w+"))==NULL)
	{
		printf("Open file sample.txt error.\n");
		return 0;
	}
	if ((fpkreg = fopen("largesamplekernelreg.txt", "w+"))==NULL)
	{
		printf("Open file samplekernelreg.txt error.\n");
		return 0;
	}
	if (LAPLACE == 0)
	{
		if ((fpkimg = fopen("largesamplekernelimg.txt", "w+"))==NULL)
		{
			printf("Open file samplekernelimg.txt error.\n");
			return 0;
		}
		if ((fpkimgreg = fopen("largesamplekernelimgreg.txt", "w+"))==NULL)
		{
			printf("Open file samplekernelimgreg.txt error.\n");
			return 0;
		}
		fprintf(fpkimg, "%d\n", size);
		fprintf(fpkimgreg, "%d\n", size);
		fprintf(fpkimg, "%.15lf %.15lf\n", samplex, sampley);
		fprintf(fpkimgreg, "%.15lf %.15lf\n", samplex, sampley);
		if (INTERIOR_MODE == 0)		//	EXTERIOR MODE, u IS OUTSIDE	//
		{
			if (cal_dist(samplex, sampley, DIST, DIST_DX, DIST_DY) <= 0.0)	//	OUTSIDE		//
			{
				fprintf(fpkimg, "%.15lf\n", cal_uimg(samplex, sampley));
				fprintf(fpkimgreg, "%.15lf\n", cal_uimg(samplex, sampley));
			}
			else					//	INSIDE		//
			{
				fprintf(fpkimg, "%.15lf\n", cal_vdlimg(samplex, sampley));
				fprintf(fpkimgreg, "%.15lf\n", cal_vdlimg(samplex, sampley));
			}
		}
		else				//	INTERIOR MODE, u IS INSIDE	//
		{
			if (cal_dist(samplex, sampley, DIST, DIST_DX, DIST_DY) >= 0.0)	//	INSIDE		//
			{
				fprintf(fpkimg, "%.15lf\n", cal_uimg(samplex, sampley));
				fprintf(fpkimgreg, "%.15lf\n", cal_uimg(samplex, sampley));
			}
			else					//	OUTSIDE		//
			{
				fprintf(fpkimg, "%.15lf\n", cal_vdlimg(samplex, sampley));
				fprintf(fpkimgreg, "%.15lf\n", cal_vdlimg(samplex, sampley));
			}
		}
	}

	
	
	
	fprintf(fpk, "%d\n", size);
	fprintf(fpkreg, "%d\n", size);
	fprintf(fpk, "%.15lf %.15lf\n", samplex, sampley);
	fprintf(fpkreg, "%.15lf %.15lf\n", samplex, sampley);
//	if (LAPLACE==0)
//	{
//		if (NONHOM_MODE==0)
//		{
//			fprintf(fpk, "%lf\n", 0.0);
//			fprintf(fpkreg, "%lf\n", 0.0);
//		}
//		else
//		{
//			cal_nonhom(samplex, sampley, tempvalue);
//			fprintf(fpk, "%.15lf\n", tempvalue[0]);
//			fprintf(fpkreg, "%.15lf\n", tempvalue[1]);
//		}
//	}

	if (INTERIOR_MODE == 0)		//	EXTERIOR MODE, u IS OUTSIDE	//
	{
		if (cal_dist(samplex, sampley, DIST, DIST_DX, DIST_DY) <= 0.0)	//	OUTSIDE		//
		{
			fprintf(fpk, "%.15lf\n", cal_u(samplex, sampley));
			fprintf(fpkreg, "%.15lf\n", cal_u(samplex, sampley));
		}
		else					//	INSIDE		//
		{
			fprintf(fpk, "%.15lf\n", cal_vdl(samplex, sampley));
			fprintf(fpkreg, "%.15lf\n", cal_vdl(samplex, sampley));
		}
	}
	else				//	INTERIOR MODE, u IS INSIDE	//
	{
		if (cal_dist(samplex, sampley, DIST, DIST_DX, DIST_DY) >= 0.0)	//	INSIDE		//
		{
			fprintf(fpk, "%.15lf\n", cal_u(samplex, sampley));
			fprintf(fpkreg, "%.15lf\n", cal_u(samplex, sampley));
		}
		else					//	OUTSIDE		//
		{
			fprintf(fpk, "%.15lf\n", cal_vdl(samplex, sampley));
			fprintf(fpkreg, "%.15lf\n", cal_vdl(samplex, sampley));
		}
	}
	
	if (ABS_MODE != 1)
		epsilon = FACTOR * H;

	for ( i = 0; i < SIZE; i++)
	{
		zy = INT_L + H * i;
		for ( j = 0; j < SIZE; j++)		
		{
			zx = INT_L + H*j;

			if (ABS_MODE == 1)
			{
				epsilon = ABS_EPSILON;
			}
			else if (ABS_MODE == 0)
			{
				thisgradx = gradientdx(zx, zy);
				thisgrady = gradientdy(zx, zy);
//				epsilon = (fabs(thisgradx) + fabs(thisgrady)) * FACTOR * H;
			}
	
//			thisdist = distance(zx, zy);
			thisdist = DIST[i][j];
			if (fabs(thisdist) > epsilon)
				continue;
			counti++;
//			zstarx = starx(zx, zy);
//			zstary = stary(zx, zy);
			zstarx = zx - thisdist * DIST_DX[i][j];
			zstary = zy - thisdist * DIST_DY[i][j];


//			gradx = cal_interpgrad(1, zstarx, zstary);
//			grady = cal_interpgrad(2, zstarx, zstary);
			cal_interpgrad(zstarx, zstary, DIST_DX, DIST_DY, &gradx, &grady);

			partialx = cal_partial(1, samplex, sampley, zstarx, zstary);
			partialy = cal_partial(2, samplex, sampley, zstarx, zstary);
			

			if (LAPLACE == 0)
			{
				partialximg = cal_partial(3, samplex, sampley, zstarx, zstary);
				partialyimg = cal_partial(4, samplex, sampley, zstarx, zstary);
			}


			//	for interior problem	//
			if (INTERIOR_MODE==1)
			{
				partial = -1.0*( gradx*partialx + grady*partialy );	//	outer normal
				if (LAPLACE == 0)
					partialimg = -1.0*( gradx*partialximg + grady*partialyimg );
			}
			
			//	for exterior problem	//
			else
			{
				partial = gradx*partialx + grady*partialy;		//	inner normal
				if (LAPLACE == 0)
					partialimg = gradx * partialximg + grady * partialyimg;
			}


			if ( ( (samplex-zstarx)*(samplex-zstarx) + (sampley-zstary)*(sampley-zstary) ) < threshold )
			{
				if (INTERIOR_MODE==1)
				{
					partial_reg = regfactor;
					if (LAPLACE == 0)
						partialimg_reg = regfactorimg;
				}
				else
				{
					partial_reg = -1.0*regfactor;
					if (LAPLACE == 0)
						partialimg_reg = -1.0*regfactorimg;
				}
			}
			else
			{
				partial_reg = partial;
				if (LAPLACE == 0)
					partialimg_reg = partialimg;
			}

			delta = cal_delta(epsilon, thisdist);
			thisJ = 1. + thisdist * CURV[i][j];
			temp = H*H* partial * delta * thisJ;
			fprintf(fpk, "%.15lf ", temp);
			temp = H*H* partial_reg * delta * thisJ;
			fprintf(fpkreg, "%.15lf ", temp);
			if (LAPLACE == 0)
			{
				temp = H*H * partialimg * delta * thisJ;
				fprintf(fpkimg, "%.15lf ", temp);
				temp = H*H * partialimg_reg * delta * thisJ;
				fprintf(fpkimgreg, "%.15lf ", temp);
			}
		}
	}
	fclose(fpk);
	fclose(fpkreg);
	if (LAPLACE == 0)
	{
		fclose(fpkimg);
		fclose(fpkimgreg);
	}
	return 0;
}

int BiCGSTAB(int size)
{
	int i, j, k, counter;
	FILE *fp3Dx, *fprealb;
	FILE *fp3DC;		//	DEBUG	//	

	double nonhomogeneous = 0.0;
	double zx, zy;
	double zstarx[size], zstary[size];
	double epsilon = FACTOR * H;
	double tempC;
	int count = 0;;

	double r[size], rhat[size], b[size], x[size], v[size], p[size], s[size], t[size];
	double error[size];
	double alpha = 1.0, rhoi = 1.0, rhoi1 = 1.0, omega = 1.0;
	double beta;
	double errornorm, preerrornorm, bnorm;
	double sum = 0.0, regsum = 0.0;
	double BC_sum = 0.0;		//	for zero BC	//
	double accpower;
	double accuracy = DEFAULTACC;

	int countrun = 0;
	int checksize;

	char denfilename[50];

//	static int flag = 0;

	time_t loop_initime, loop_endtime, seconds;

	sprintf(denfilename, "DensityN%d.txt", (int)(N));
		
	if ((fp3Dx = fopen( denfilename, "w+"))==NULL)
	{
		printf("Open file %s error.\n", denfilename);
		return 0;
	}
//	fprintf(fp3Dx, "%d\n", size);

	if ((fprealb = fopen("largerealb.txt", "r"))==NULL)
	{
		printf("Open file large3Drealb.txt error.\n");
		return 0;
	}
	if (TEST_CONDITION==1)
	{
		if ((fp3DC = fopen("large3DC.txt", "w+"))==NULL)
		{
			printf("Open file large3DC.txt error.\n");
			return 0;
		}
	}
	//	fprintf(fp3DC, "%d\n", size);

	if (INTERACTIVE_ACCURACY == 1) 
	{
		printf("Enter the accuracy in 10^x ( 0 > x > -10)\n");
		scanf("%lf", &accpower);
		if ( (accpower > -10.0) && (accpower < 0.0) )
			accuracy = pow(10, accpower);
	}
	printf("Using accuracy = %.8E.\t", accuracy);

	//	INITIALIZING THE THINGS NEEDED FOR BICONJUGATE GRADIENT STABILIZED	//

	fscanf(fprealb, "%d", &checksize);
	if (checksize != size)
	{
		printf("size error.\n");
		return 0;
	}
	for ( counter = 0; counter < size; counter++)
	{
		fscanf(fprealb, "%lf", &b[counter]);
		if (b[counter]!=b[counter])
			printf("b[%d] = %lf.\n", counter, b[counter]);
		r[counter] = b[counter];
		v[counter] = 0.0;
		p[counter] = 0.0;
		rhat[counter] = r[counter];
		if (rhat[counter] != r[counter])
			printf("counter = %d.\n, rhat = %lf, r = %lf.\n", counter, rhat[counter], r[counter]);
	}
	fclose(fprealb);
	bnorm = norm(size, b);
	///////////////////////////////////////////////////////////////////////////////////////

//	if (NONHOM_MODE==1)
///	{
//		nonhomogeneous = cal_nonhom();
//		count = 0;
//		for (int i = 0; i < N; i++)
//		{
//			zy = INT_L + h * double(i);
//			for (int j = 0; j < N; j++)
//			{
//				zx = INT_L + h * double (j);
//				if ( fabs(distance( zx, zy)) < epsilon )
///				{
//					zstarx[count] = starx(zx, zy);
///					zstary[count] = stary(zx, zy);
//					count++;
//				}
//			}
//		}
//
//		for (int counter = 0; counter < size; counter++)
//		{
//			b[counter] -= cal_nonhom(zstarx[counter], zstary[counter]);
//			r[counter] = b[counter];
//			rhat[counter] = r[counter];
//			if ( (b[counter]!=b[counter]) || (b[counter]+1.0==b[counter]) )
//			{
//				printf("b[%d] is %lf.\n", counter, b[counter]);
//				exit(0);
//			}
//		}
//		
//	}

	for ( i = 0; i < size; i++)
	{
		BC_sum += fabs(r[i]);
	}
	if ( BC_sum < pow(10,-13.0))	//	for zero BC	//
	{
		printf("Here.\n");
		r[0] = 0.000000001;
		rhat[0] = 0.000000001;
	}

	if ( BC_sum < pow(10, -15.0) )
	{
		printf("Zero BC.\n");
		for ( i = 0; i < size; i++)
		{
			fprintf(fp3Dx, "%.15lf ", 0.0);
			
		}
		fclose(fp3Dx);
		return 0;
	}


	for ( i = 1; i <= 20000; i++)
	{
		loop_initime = time(NULL);

		rhoi1 = rhoi;					//	

		rhoi = inner(size, rhat, r);			//	rho = (r^, r)			//
//		printf("rhoi = %lf, rhat[0] = %lf, r[0] = %lf.\n", rhoi, rhat[0], r[0]);
//		if (rhoi!=rhoi)
//			printf("rhoi.\n");
		beta = (rhoi*alpha)/(rhoi1*omega);		//	beta = (rhoi/rhoi1) * (a/w)	//
		//		printf("beta = %lf\n", beta);
		for ( j = 0; j < size; j++)
		{
			p[j] = r[j] + beta * (p[j] - omega * v[j]);	//	p = r + beta( p - w * v )  //
//			if ( (p[j] != p[j])|| (p[j] + 1.0 ==p[j]) )
//				printf("j = %d, r[j] = %lf, v[j] = %lf, beta = %lf.\n", j, r[j], v[j], beta);
		}

		//////	FIRST MATRIX MULTIPLICATION	////////////
		for ( j = 0; j < size; j++)
		{
			sum = 0.0;
			for ( k = 0; k < size; k++)
			{
				tempC = get_C(size, j, k);
				if (i==1)
				{
					if (TEST_CONDITION==1)
						fprintf(fp3DC, "%.8lf ", tempC);
				}
				sum += (tempC * p[k]);
			}
			if (i==1)
			{
				if (TEST_CONDITION==1)
					fprintf(fp3DC, "\n");
			}
			v[j] = sum;				//	v = A * p	//
			if ((v[j] != v[j])||(v[j] +1 == v[j]))
			{
				printf("v wrong, j = %d, v[j] = %lf.\n", j, v[j]);
				exit(0);
			}
//			printf("v[%d] = %.15lf\n", j, v[j]);
		}
		loop_endtime = time(NULL);
		seconds = loop_endtime - loop_initime;
//		printf("Made to 1st. Took %ld minutes %ld seconds.\n", seconds/60, seconds % 60);
		/////	FIRST MATRIX MULTIPLICATION ENDS	////


		if (TEST_CONDITION==1)
		{
			fclose(fp3DC);
			printf("Printed matrix for condition number.\n");
			exit(0);
		}



//		printf("v[0] = %.15lf.\n", v[0]);
//		printf("rhoi = %.15lf.\n", rhoi);
		alpha = rhoi/inner(size, rhat, v);	//	a = rho/ (r^, v)	//
		for ( j = 0; j < size; j++)		
		{
			s[j] = r[j] - alpha * v[j];	//	s = r - a * v	//
			if ((s[j] !=s[j])||(s[j] + 1 ==s[j]))
			{
				printf("s wrong, j = %d, alpha = %lf, r[j] = %lf\n", j, alpha, r[j]);
				exit(0);
			}
		}

		//////	SECOND MATRIX MULTIPLICATION	/////////////
		for ( j = 0; j < size; j++)
		{
			sum = 0.0;
			for ( k = 0; k < size; k++)
			{
				sum += (get_C(size, j, k) * s[k]);
			}
			t[j] = sum;				//	t = A * s		//
			if (t[j]!=t[j])
			{
				printf("t wrong, j = %d.\n", j);
				printf("alpha = %lf, r[%d] = %lf.\n", alpha, j, r[j]);
				exit(0);
			}
		}
		///////	SECOND MATRIX MULTIPLICATION ENDS	//////

		omega = inner(size, t, s)/inner(size, t, t);	//	w = (t, s)/(t, t)	//
		for ( j = 0; j < size; j++)
		{
			x[j] = x[j] + alpha * p[j] + omega * s[j];	//	x = x + a * p + w * s	//
		}
	
		///////	THIRD MATRIX MULTIPLICATION	/////////////
//		for (int j = 0; j < size; j++)
//		{
//			sum = 0.0;
//			for (int k = 0; k < size; k++)
//			{
//				sum += (get_C(size, j, k) * x[k]);
//			}
//			error[j] = sum - b[j];			//	err = b - A * x		//
//			if (error[j] != error[j])
//				printf("error wrong, j = %d, x = %lf, b = %lf, omega = %lf.", j, x[j], b[j], omega);
//		}

		for ( j = 0; j < size; j++)
		{
			r[j] = s[j] - omega * t[j];		//	r = s - w * t		//
		}

		preerrornorm = errornorm;
//		errornorm = norm(size, error);
		errornorm = norm(size,r)/bnorm;
		printf("Run %d, error = %.8E\n", i, errornorm);


		loop_endtime = time(NULL);
		seconds = loop_endtime - loop_initime;
//		printf("This run took %ld minutes and %ld seconds.\n", seconds/60, seconds % 60);
		countrun++;

		if (errornorm < accuracy)
			break;
		if (errornorm != errornorm)
			break;
		if (fabs(preerrornorm - errornorm) < pow(10, -3) * accuracy)
			break;
		

	}

	printf("After %d runs, the error is %.10E.\n", countrun, errornorm);
	
	for ( i = 0; i < size; i++)
	{
		fprintf(fp3Dx, "%.15lf ", x[i]);
		Density[i] = x[i];
	}
	
	fclose(fp3Dx);
//	fclose(fp3DC);
	return 0;
}


int BiCGcomplex(int size, double epsilon)
{
	int i, j, k, counter;
	FILE *fp3Dx, *fprealb, *fp3Dximg, *fpimgb;
	FILE *fpverifyb, *fpverifybimg;
	FILE *fp3DC;		//	DEBUG	//	
	FILE *fp3DCimg;		//	DEBUG	//

	int compsize = 2 * size;
	double nonhomogeneous = 0.0;
	double nonhomogeneousimg = 0.0;
	double zx, zy;
	double zstarx[size], zstary[size];
//	double epsilon = FACTOR * h;
	int count = 0;;

	double tempC[2];
//	double r[size], rhat[size], x[size], v[size], p[size], s[size], t[size];
//	double rimg[size], rhatimg[size], ximg[size], vimg[size], pimg[size], simg[size], timg[size];
	double x[size], ximg[size], xtilde[size], xtildeimg[size], r[size], rimg[size], rtilde[size], rtildeimg[size];
	double p[size], pimg[size], ptilde[size], ptildeimg[size], pold[size], poldimg[size], ptildeold[size], ptildeoldimg[size];
	double Akp[size], Akpimg[size], Akptilde[size], Akptildeimg[size];

//	double b[3] = TESTB;
//	double bimg[3] = TESTBIMG;
	double b[size], bimg[size];




	double error[size], errorimg[size];
	double alpha = 1.0, rhoi = 1.0, rhoi1 = 1.0, omega = 1.0;
	double alphaimg = 0.0, rhoiimg = 0.0, rhoi1img = 0.0, omegaimg = 0.0;
	double beta, betaimg;
	double errornorm, preerrornorm;
	double comptemp, comptempnum, comptempnumimg, comptempden, comptempa, comptempb, comptempc, comptempd;
	double innerrold, innerroldimg;

	double sum = 0.0, regsum = 0.0, sumimg = 0.0, regsumimg = 0.0;
	double normb = 0.0;
	double BC_sum = 0.0;		//	for zero BC	//
	double accpower;
	double accuracy = DEFAULTACC;

	int countrun = 0;
	int checksize, checksizeimg;

	char denfilename[50], denimgfilename[50];

//	static int flag = 0;

	time_t loop_initime, loop_endtime, seconds;

	sprintf(denfilename, "DensityN%dk%d.txt", (int)(N), (int)(WAVE));
	if ((fp3Dx = fopen( denfilename, "w+"))==NULL)
	{
		printf("Open file %s error.\n", denfilename);
		return 0;
	}
	sprintf(denimgfilename, "DensityimgN%dk%d.txt", (int)(N), (int)(WAVE));
	if ((fp3Dximg = fopen( denimgfilename, "w+"))==NULL)
	{
		printf("Open file %s error.\n", denimgfilename);
		return 0;
	}
//	fprintf(fp3Dx, "%d\n", size);

	if ((fpverifyb = fopen("verifyb.txt", "w+"))==NULL)
	{
		printf("Open file verifyb.txt error.\n");
		exit(0);
	}
	if ((fpverifybimg = fopen("verifybimg.txt", "w+"))==NULL)
	{
		printf("Open file verifybimg.txt error.\n");
		exit(0);
	}
	
	if (TEST_CONDITION==1)
	{
		if ((fp3DC = fopen("testC.txt", "w+"))==NULL)
		{
			printf("Open file testC.txt error.\n");
			return 0;
		}
	
//
//
//	fprintf(fp3DC, "%d\n", size);
	
		if ((fp3DCimg = fopen("testCimg.txt", "w+"))==NULL)
		{
			printf("Open file testCimg.txt error.\n");
			return 0;
		}
//	fprintf(fp3DCimg, "%d\n", size);
	}

	if (INTERACTIVE_ACCURACY == 1)
	{
		printf("Enter the accuracy in 10^x ( 0 > x > -10)\n");
		scanf("%lf", &accpower);
		if ( (accpower > -10.0) && (accpower < 0.0) )
			accuracy = pow(10, accpower);
		printf("Using accuracy = %.8E.\n", accuracy);
	}


	//	INITIALIZING THE THINGS NEEDED FOR BICONJUGATE GRADIENT STABILIZED	//


	if ((fprealb = fopen("largerealb.txt", "r"))==NULL)
	{
		printf("Open file largerealb.txt error.\n");
		return 0;
	}
	if ((fpimgb = fopen("largeimgb.txt", "r"))==NULL)
	{
		printf("Open file largeimgb.txt error.\n");
		return 0;
	}

	fscanf(fprealb, "%d", &checksize);
	fscanf(fpimgb, "%d", &checksizeimg);
	if ((checksize != size) || (checksizeimg!=size))
	{
		printf("size error. Realsize = %d, realpart = %d, imgpart = %d\n", size, checksize, checksizeimg);
		return 0;
	}
	for ( counter = 0; counter < size; counter++)
	{
		fscanf(fprealb, "%lf", &b[counter]);
		fscanf(fpimgb, "%lf", &bimg[counter]);

		if ((b[counter]!=b[counter])||(b[size+counter]!=b[size+counter]))
		{
			printf("b[%d] = %lf.\n", counter, b[counter]);
			printf("bimg[%d] = %lf.\n", counter, bimg[counter]);
			exit(0);
		}
	}
	fclose(fprealb);
	fclose(fpimgb);


////////////////////////////////////////////////////////////////////////////
//	TESTB;	
//	bimg[3] = TESTBIMG;

	fprintf(fpverifyb, "%d\n", size);
	fprintf(fpverifybimg, "%d\n", size);
	for (i = 0; i < size ; i++)
	{
		fprintf(fpverifyb, "%lf ", b[i]);
		fprintf(fpverifybimg, "%lf ", bimg[i]);
	}
////////////////////////////////////////////////////////////////////////////


	normb = sqrt(inner(size, b, b) + inner(size, bimg, bimg));
	printf("b norm is %lf\t", normb);

	for ( counter = 0; counter < size; counter++)
	{
		x[counter] = 0.0;
		ximg[counter] = 0.0;
		r[counter] = b[counter];
		rimg[counter] = bimg[counter];
		p[counter] = r[counter];
		pimg[counter] = rimg[counter];

		rtilde[counter] = r[counter];
		rtildeimg[counter] = (-1.0)*rimg[counter];
		ptilde[counter] = rtilde[counter];
		ptildeimg[counter] = rtildeimg[counter];

		if ( (rtilde[counter] != rtilde[counter]) || (rtildeimg[counter] != rtildeimg[counter]) )
		{
			printf("counter = %d.\n, rtilde = %lf, r = %lf.\n", counter, rtilde[counter], r[counter]);
			printf("rtildeimg = %lf, rimg = %lf.\n", rtildeimg[counter], rimg[counter]);
		}
	}
	//	(r, r^)	//
	innerrold = inner(size, r, rtilde) + inner(size, rimg, rtildeimg);
	innerroldimg = inner(size, rimg, rtilde) - inner(size, r, rtildeimg);

	///////////////////////////////////////////////////////////////////////////////////////


	for ( i = 0; i < size; i++)
	{
		BC_sum += fabs(r[i]);
		BC_sum += fabs(rimg[i]);
	}
	if ( BC_sum < pow(10,-13.0))	//	for zero BC	//
	{
		printf("Here.\n");
		r[0] = 0.000000001;
		rimg[0] = 0.000000001;
		rtilde[0] = 0.000000001;
		rtildeimg[0] = 0.000000001;
	}

	if ( BC_sum < pow(10, -15.0) )
	{
		printf("Zero BC.\n");
		for ( i = 0; i < size; i++)
		{
			fprintf(fp3Dx, "%.15lf ", 0.0);
			fprintf(fp3Dximg, "%.15lf ", 0.0);
		}
		fclose(fp3Dx);
		fclose(fp3Dximg);
		return 0;
	}


	for ( i = 1; i <= 200; i++)
	{
		loop_initime = time(NULL);

		//////	FIRST MATRIX MULTIPLICATION	////////////
		for ( j = 0; j < size; j++)
		{
			sum = 0.0;
			sumimg = 0.0;
			for ( k = 0; k < size; k++)
			{
				get_Ccomp(size, j, k, epsilon, tempC);
				if (i==1)
				{
					if (TEST_CONDITION==1)
					{
						fprintf(fp3DC, "%lf ", tempC[0]);
						fprintf(fp3DCimg, "%lf ", tempC[1]);
					}
				}
				sum += (tempC[0] * p[k] - tempC[1] * pimg[k]);		//	A * p	//
				sumimg += (tempC[0] * pimg[k] + tempC[1] * p[k]);	//	A * p	//
			}
			if (i==1)
			{
				if (TEST_CONDITION==1)
				{
					fprintf(fp3DC, "\n");
					fprintf(fp3DCimg, "\n");
				}
			}
			Akp[j] = sum;				//	v = A * p	//
			Akpimg[j] = sumimg;

//			if ((Akp[j] != Akp[j])||(Akp[j] +1 == Akp[j])||(Akpimg[j] != Akpimg[j])||(Akpimg[j]+1==Akpimg[j]))
//			{
//				printf("Akp wrong, j = %d, Akp[j] = %lf, Akpimg[j] = %lf.\n", j, Akp[j], Akpimg[j]);
//				exit(0);
//			}
//			printf("v[%d] = %.15lf\n", j, v[j]);
		}
		loop_endtime = time(NULL);
		seconds = loop_endtime - loop_initime;

//		if (i == 1)
//			printf("Made to 1st. Took %ld minutes %ld seconds.\n", seconds/60, seconds % 60);
		/////	FIRST MATRIX MULTIPLICATION ENDS	////

		if (TEST_CONDITION==1)
		{
			fclose(fp3DC);
			fclose(fp3DCimg);
			printf("Printed matrices for condition number test.\n");
			exit(0);
		}





		//////	SECOND MATRIX MULTIPLICATION	/////////////
		for ( j = 0; j < size; j++)
		{
			sum = 0.0;
			sumimg = 0.0;
			for ( k = 0; k < size; k++)
			{
				get_Ccomp(size, k, j, epsilon, tempC);
				sum += (tempC[0]*ptilde[k] + tempC[1]*ptildeimg[k]);
				sumimg += (tempC[0]*ptildeimg[k] - tempC[1]*ptilde[k]);
			}
			Akptilde[j] = sum;				//	Akptilde = A^H * p		//
			Akptildeimg[j] = sumimg;
//			if ((Akptilde[j]!=Akptilde[j])||(Akptildeimg[j]!=Akptildeimg[j]))
//			{
//				printf("Akptilde wrong, j = %d, Akptilde[j] = %lf, Akptildeimg[j] = %lf.\n", \
//				j, Akptilde[j], Akptildeimg[j]);
//				printf("alpha = %lf, r[%d] = %lf.\n", alpha, j, r[j]);
//				exit(0);
//			}
		}
		///////	SECOND MATRIX MULTIPLICATION ENDS	//////


		//	DEFINE IF A = a+bi, B = c+di, then (A,B) = sigma(A*B') = (a+bi)dot(c-di)	//
		//	THAT IS, CONJUGATE THE LATTER VECTOR	//

		comptempa = inner(size, r,rtilde) + inner(size, rimg, rtildeimg);
		comptempb = inner(size, rimg,rtilde) - inner(size, r, rtildeimg);
		comptempc = inner(size, Akp, ptilde) + inner(size, Akpimg, ptildeimg);
		comptempd = inner(size, Akpimg, ptilde) - inner(size, Akp, ptildeimg);
		comptempnum = comptempa*comptempc + comptempb*comptempd;
		comptempnumimg = comptempb*comptempc - comptempa*comptempd;
		comptempden = comptempc*comptempc + comptempd*comptempd;
		alpha = comptempnum/comptempden;
		alphaimg = comptempnumimg/comptempden;

		for ( j = 0; j < size; j++)
		{
			x[j] = x[j] + alpha * p[j] - alphaimg * pimg[j];
			ximg[j] = ximg[j] + alphaimg * p[j] + alpha * pimg[j];
			r[j] = r[j] - alpha * Akp[j] + alphaimg * Akpimg[j];
			rimg[j] = rimg[j] - alphaimg * Akp[j] - alpha * Akpimg[j];

			rtilde[j] = rtilde[j] - alpha * Akptilde[j] - alphaimg * Akptildeimg[j];
			rtildeimg[j] = rtildeimg[j] + alphaimg * Akptilde[j] - alpha * Akptildeimg[j];
			
//			if ( (p[j] != p[j])|| (p[j] + 1.0 ==p[j]) || (alpha != alpha)  || (alphaimg != alphaimg) || \
//				(pimg[j] != pimg[j]) )
//			{
//				printf("j = %d, r = (%lf, %lf), p = (%lf, %lf),alpha = (%lf, %lf).\n", \
//				j, r[j], rimg[j], p[j], pimg[j], alpha, alphaimg);
//				exit(0);
//			}
		}


		//	BEGIN beta = (r^,r)_k+1/(r^, r)k		//
		comptempc = innerrold;
		comptempd = innerroldimg;	
		comptempa = inner(size, r, rtilde) + inner(size, rimg, rtildeimg);
		comptempb = inner(size, rimg, rtilde) - inner(size, r, rtildeimg);
		comptempnum = comptempa*comptempc + comptempb*comptempd;
		comptempnumimg = comptempb*comptempc - comptempa*comptempd;
		comptempden = comptempc*comptempc + comptempd*comptempd;
		beta = comptempnum/comptempden;
		betaimg = comptempnumimg/comptempden;

		innerrold = comptempa;
		innerroldimg = comptempb;

		//	END a = rho/ (r^, v)	//

		for ( j = 0; j < size; j++)		
		{
			pold[j] = p[j];
			poldimg[j] = pimg[j];
			ptildeold[j] = ptilde[j];
			ptildeoldimg[j] = ptildeimg[j];
		}

		for ( j = 0; j < size; j++)
		{
			p[j] = r[j] + beta * pold[j] - betaimg * poldimg[j];		//	p = r + beta * p	//
			pimg[j] = rimg[j] + betaimg * pold[j] + beta * poldimg[j];

			ptilde[j] = rtilde[j] + beta * ptildeold[j] + betaimg * ptildeoldimg[j];
			ptildeimg[j] = rtildeimg[j] - betaimg * ptildeold[j] + beta * ptildeoldimg[j];


			if ((p[j] !=p[j])||(p[j] + 1 ==p[j])||(pimg[j] !=pimg[j])||(pimg[j]+1==pimg[j]))
			{
				printf("s wrong, j = %d, alpha = %lf, r[j] = %lf, rimg[j] = %lf.\n", j, alpha, r[j], rimg[j]);
				exit(0);
			}
		}



		preerrornorm = errornorm;
		errornorm = sqrt(inner(size, r, r) + inner(size, rimg, rimg))/normb;
//		printf("Run %d, error = %.8E\n", i, errornorm);


		loop_endtime = time(NULL);
		seconds = loop_endtime - loop_initime;
//		if (i == 1)
//			printf("This run took %ld minutes and %ld seconds.\n", seconds/60, seconds % 60);
		countrun++;

		if (errornorm < accuracy)
			break;
		if (errornorm != errornorm)
			break;
		if (fabs(preerrornorm - errornorm) < pow(10, -3) * accuracy)
			break;
		
	}

	printf("After %d runs, the error is %.10E.\n", countrun, errornorm);
	
	for ( i = 0; i < size; i++)
	{
		fprintf(fp3Dx, "%.15lf ", x[i]);
		fprintf(fp3Dximg, "%.15lf ", ximg[i]);
		if ( (x[i] != x[i])|| (ximg[i] != ximg[i]) )
		{
			printf("AT i = %d, x = (%lf, %lf)\n", i, x[i], ximg[i]);
			exit(0);
		}

		Density[i] = x[i];
		Densityimg[i] = ximg[i];
	}
	get_Ccomp(size, -1, -1, epsilon, tempC);
	fclose(fp3Dx);
	fclose(fp3Dximg);
//	fclose(fp3DC);
//	fclose(fp3DCimg);
	return 0;
}


int parallel_BiCGcomplex(int size, double epsilon, int chunk)
{
	int i, j, k, counter;
	FILE *fp3Dx, *fprealb, *fp3Dximg, *fpimgb;
	FILE *fpverifyb, *fpverifybimg;
	FILE *fp3DC;		//	DEBUG	//	
	FILE *fp3DCimg;		//	DEBUG	//

	int compsize = 2 * size;
	double nonhomogeneous = 0.0;
	double nonhomogeneousimg = 0.0;
//	double zx, zy;
//	double zstarx[size], zstary[size];
	int count = 0;;

	double tempC[2];
	double *x, *ximg, *r, *rimg, *rtilde, *rtildeimg;
	double *p, *pimg, *ptilde, *ptildeimg, *pold, *poldimg, *ptildeold, *ptildeoldimg;
	double *Akp, *Akpimg, *Akptilde, *Akptildeimg;

//	double b[3] = TESTB;
//	double bimg[3] = TESTBIMG;
	double *b, *bimg;
	double **matrix, **matriximg;



	double *error, *errorimg;
	double alpha = 1.0, rhoi = 1.0, rhoi1 = 1.0, omega = 1.0;
	double alphaimg = 0.0, rhoiimg = 0.0, rhoi1img = 0.0, omegaimg = 0.0;
	double beta, betaimg;
	double errornorm, preerrornorm;
	double comptemp, comptempnum, comptempnumimg, comptempden, comptempa, comptempb, comptempc, comptempd;
	double innerrold, innerroldimg;

	double sum = 0.0, regsum = 0.0, sumimg = 0.0, regsumimg = 0.0;
	double normb = 0.0;
	double BC_sum = 0.0;		//	for zero BC	//
	double accpower;
	double accuracy = DEFAULTACC;

	int countrun = 0;
	int checksize, checksizeimg;

	char denfilename[50], denimgfilename[50];

//	static int flag = 0;

	time_t loop_initime, loop_endtime, seconds;


	x = (double *) malloc(size * sizeof(double));
	r = (double *) malloc(size * sizeof(double));
	p = (double *) malloc(size * sizeof(double));
	pold = (double *) malloc(size * sizeof(double));
	rtilde = (double *) malloc(size * sizeof(double));
	ptilde = (double *) malloc(size * sizeof(double));
	ptildeold = (double *) malloc(size * sizeof(double));
	ximg = (double *) malloc(size * sizeof(double));
	rimg = (double *) malloc(size * sizeof(double));
	pimg = (double *) malloc(size * sizeof(double));
	poldimg = (double *) malloc(size * sizeof(double));
	rtildeimg = (double *) malloc(size * sizeof(double));
	ptildeimg = (double *) malloc(size * sizeof(double));
	ptildeoldimg = (double *) malloc(size * sizeof(double));
	b = (double *) malloc(size * sizeof(double));
	bimg = (double *) malloc(size * sizeof(double));
	error = (double *) malloc(size * sizeof(double));
	errorimg = (double *) malloc(size * sizeof(double));
	Akp = (double *) malloc(size * sizeof(double));
	Akptilde = (double *) malloc(size * sizeof(double));
	Akpimg = (double *) malloc(size * sizeof(double));
	Akptildeimg = (double *) malloc(size * sizeof(double));



	matrix = (double **) malloc(size * sizeof(double *));
	matriximg = (double **) malloc(size * sizeof(double *));
	for (i = 0; i < size; i++)
	{
		matrix[i] = (double *) malloc(size * sizeof(double));
		matriximg[i] = (double *) malloc(size * sizeof(double));
	}




	sprintf(denfilename, "DensityN%dk%d.txt", (int)(N), (int)(WAVE));
	if ((fp3Dx = fopen( denfilename, "w+"))==NULL)
	{
		printf("Open file %s error.\n", denfilename);
		return 0;
	}
	sprintf(denimgfilename, "DensityimgN%dk%d.txt", (int)(N), (int)(WAVE));
	if ((fp3Dximg = fopen( denimgfilename, "w+"))==NULL)
	{
		printf("Open file %s error.\n", denimgfilename);
		return 0;
	}
//	fprintf(fp3Dx, "%d\n", size);

	if ((fpverifyb = fopen("verifyb.txt", "w+"))==NULL)
	{
		printf("Open file verifyb.txt error.\n");
		exit(0);
	}
	if ((fpverifybimg = fopen("verifybimg.txt", "w+"))==NULL)
	{
		printf("Open file verifybimg.txt error.\n");
		exit(0);
	}
	
	if (TEST_CONDITION==1)
	{
		if ((fp3DC = fopen("testC.txt", "w+"))==NULL)
		{
			printf("Open file testC.txt error.\n");
			return 0;
		}
	
//
//
//	fprintf(fp3DC, "%d\n", size);
	
		if ((fp3DCimg = fopen("testCimg.txt", "w+"))==NULL)
		{
			printf("Open file testCimg.txt error.\n");
			return 0;
		}
//	fprintf(fp3DCimg, "%d\n", size);
	}

	if (INTERACTIVE_ACCURACY == 1)
	{
		printf("Enter the accuracy in 10^x ( 0 > x > -10)\n");
		scanf("%lf", &accpower);
		if ( (accpower > -10.0) && (accpower < 0.0) )
			accuracy = pow(10, accpower);
		printf("Using accuracy = %.8E.\n", accuracy);
	}


	//	INITIALIZING THE THINGS NEEDED FOR BICONJUGATE GRADIENT STABILIZED	//


	if ((fprealb = fopen("largerealb.txt", "r"))==NULL)
	{
		printf("Open file largerealb.txt error.\n");
		return 0;
	}
	if ((fpimgb = fopen("largeimgb.txt", "r"))==NULL)
	{
		printf("Open file largeimgb.txt error.\n");
		return 0;
	}

	fscanf(fprealb, "%d", &checksize);
	fscanf(fpimgb, "%d", &checksizeimg);
	if ((checksize != size) || (checksizeimg!=size))
	{
		printf("size error. Realsize = %d, realpart = %d, imgpart = %d\n", size, checksize, checksizeimg);
		return 0;
	}
	for ( counter = 0; counter < size; counter++)
	{
		fscanf(fprealb, "%lf", &b[counter]);
		fscanf(fpimgb, "%lf", &bimg[counter]);

		if ((b[counter]!=b[counter])||(b[size+counter]!=b[size+counter]))
		{
			printf("b[%d] = %lf.\n", counter, b[counter]);
			printf("bimg[%d] = %lf.\n", counter, bimg[counter]);
			exit(0);
		}
	}
	fclose(fprealb);
	fclose(fpimgb);


////////////////////////////////////////////////////////////////////////////
//	TESTB;	
//	bimg[3] = TESTBIMG;

	fprintf(fpverifyb, "%d\n", size);
	fprintf(fpverifybimg, "%d\n", size);
	for (i = 0; i < size ; i++)
	{
		fprintf(fpverifyb, "%lf ", b[i]);
		fprintf(fpverifybimg, "%lf ", bimg[i]);
	}
////////////////////////////////////////////////////////////////////////////


	normb = sqrt(omp_inner(size, b, b, chunk) + omp_inner(size, bimg, bimg, chunk));
	printf("b norm is %lf\t", normb);

	for ( counter = 0; counter < size; counter++)
	{
		x[counter] = 0.0;
		ximg[counter] = 0.0;
		r[counter] = b[counter];
		rimg[counter] = bimg[counter];
		p[counter] = r[counter];
		pimg[counter] = rimg[counter];

		rtilde[counter] = r[counter];
		rtildeimg[counter] = (-1.0)*rimg[counter];
		ptilde[counter] = rtilde[counter];
		ptildeimg[counter] = rtildeimg[counter];

		if ( (rtilde[counter] != rtilde[counter]) || (rtildeimg[counter] != rtildeimg[counter]) )
		{
			printf("counter = %d.\n, rtilde = %lf, r = %lf.\n", counter, rtilde[counter], r[counter]);
			printf("rtildeimg = %lf, rimg = %lf.\n", rtildeimg[counter], rimg[counter]);
		}
	}
	//	(r, r^)	//
	innerrold = omp_inner(size, r, rtilde, chunk) + omp_inner(size, rimg, rtildeimg, chunk);
	innerroldimg = omp_inner(size, rimg, rtilde, chunk) - omp_inner(size, r, rtildeimg, chunk);

	///////////////////////////////////////////////////////////////////////////////////////


	for ( i = 0; i < size; i++)
	{
		BC_sum += fabs(r[i]);
		BC_sum += fabs(rimg[i]);
	}
	if ( BC_sum < pow(10,-13.0))	//	for zero BC	//
	{
		printf("Here.\n");
		r[0] = 0.000000001;
		rimg[0] = 0.000000001;
		rtilde[0] = 0.000000001;
		rtildeimg[0] = 0.000000001;
	}

	if ( BC_sum < pow(10, -15.0) )
	{
		printf("Zero BC.\n");
		for ( i = 0; i < size; i++)
		{
			fprintf(fp3Dx, "%.15lf ", 0.0);
			fprintf(fp3Dximg, "%.15lf ", 0.0);
		}
		fclose(fp3Dx);
		fclose(fp3Dximg);
		return 0;
	}


	for ( i = 1; i <= 200; i++)
	{
		loop_initime = time(NULL);

		#pragma omp parallel for private(j) schedule(static, chunk)
		for (j = 0; j < size; j++)
		{
			Akp[j] = 0.0;
			Akpimg[j] = 0.0;
			Akptilde[j] = 0.0;
			Akptildeimg[j] = 0.0;
		}

		//	FIRST OF FIRST MATRIX MULTIPLICATION	//
		if (i == 1)
		{
			for (j = 0; j < size; j++)
			{
			for (k = 0; k < size; k++)
			{
				omp_get_Ccomp(size, j ,k, epsilon, tempC);
				matrix[j][k] = tempC[0];
				matriximg[j][k] = tempC[1];
				Akp[j] += (tempC[0] * p[k] - tempC[1] * pimg[k]);		//	A * p	//
				Akpimg[j] += (tempC[0] * pimg[k] + tempC[1] * p[k]);	//	A * p	//
				if (TEST_CONDITION==1)
				{
					fprintf(fp3DC, "%lf ", tempC[0]);
					fprintf(fp3DCimg, "%lf ", tempC[1]);
				}
			}
				if (TEST_CONDITION==1)
				{
					fprintf(fp3DC, "\n");
					fprintf(fp3DCimg, "\n");
				}
			}
//			loop_endtime = time(NULL);
//			seconds = loop_endtime - loop_initime;
//			printf("Made to 1st. Took %ld minutes %ld seconds.\n", seconds/60, seconds % 60);
			if (TEST_CONDITION==1)
			{
				fclose(fp3DC);
				fclose(fp3DCimg);
				printf("Printed matrices for condition number test.\n");
				exit(0);
			}
			omp_get_Ccomp(size, -1, -1, epsilon, tempC);
		}
		else
		{
			//////	FIRST MATRIX MULTIPLICATION	////////////
			#pragma omp parallel for private(j, k) schedule(static, chunk)
			for ( j = 0; j < size; j++)
			{
			for ( k = 0; k < size; k++)
			{
				Akp[j] += (matrix[j][k] * p[k] - matriximg[j][k] * pimg[k]);		//	A * p	//
				Akpimg[j] += (matrix[j][k] * pimg[k] + matriximg[j][k] * p[k]);	//	A * p	//
			}

//			if ((Akp[j] != Akp[j])||(Akp[j] +1 == Akp[j])||(Akpimg[j] != Akpimg[j])||(Akpimg[j]+1==Akpimg[j]))
//			{
//				printf("Akp wrong, j = %d, Akp[j] = %lf, Akpimg[j] = %lf.\n", j, Akp[j], Akpimg[j]);
//				exit(0);
//			}
//			printf("v[%d] = %.15lf\n", j, v[j]);
			}
		}
		/////	FIRST MATRIX MULTIPLICATION ENDS	////


		//////	SECOND MATRIX MULTIPLICATION	/////////////
		#pragma omp parallel for private(j, k) schedule(static, chunk)
		for ( j = 0; j < size; j++)
		{
		for ( k = 0; k < size; k++)
		{
				Akptilde[j] += (matrix[k][j]*ptilde[k] + matriximg[k][j]*ptildeimg[k]);
				Akptildeimg[j] += (matrix[k][j]*ptildeimg[k] - matriximg[k][j]*ptilde[k]);
		}
//		if ((Akptilde[j]!=Akptilde[j])||(Akptildeimg[j]!=Akptildeimg[j]))
//		{
//			printf("Akptilde wrong, j = %d, Akptilde[j] = %lf, Akptildeimg[j] = %lf.\n", \
//			j, Akptilde[j], Akptildeimg[j]);
//			printf("alpha = %lf, r[%d] = %lf.\n", alpha, j, r[j]);
//			exit(0);
//		}
		}
		///////	SECOND MATRIX MULTIPLICATION ENDS	//////


		//	DEFINE IF A = a+bi, B = c+di, then (A,B) = sigma(A*B') = (a+bi)dot(c-di)	//
		//	THAT IS, CONJUGATE THE LATTER VECTOR	//

		comptempa = omp_inner(size, r,rtilde, chunk) + omp_inner(size, rimg, rtildeimg, chunk);
		comptempb = omp_inner(size, rimg,rtilde, chunk) - omp_inner(size, r, rtildeimg, chunk);
		comptempc = omp_inner(size, Akp, ptilde, chunk) + omp_inner(size, Akpimg, ptildeimg, chunk);
		comptempd = omp_inner(size, Akpimg, ptilde, chunk) - omp_inner(size, Akp, ptildeimg, chunk);
		comptempnum = comptempa*comptempc + comptempb*comptempd;
		comptempnumimg = comptempb*comptempc - comptempa*comptempd;
		comptempden = comptempc*comptempc + comptempd*comptempd;
		alpha = comptempnum/comptempden;
		alphaimg = comptempnumimg/comptempden;

		#pragma omp parallel sections private(j)
		{
			#pragma omp section
			{
				for ( j = 0; j < size; j++)
					x[j] = x[j] + alpha * p[j] - alphaimg * pimg[j];
			}
			#pragma omp section
			{
				for ( j = 0; j < size; j++)
					ximg[j] = ximg[j] + alphaimg * p[j] + alpha * pimg[j];
			}
			#pragma omp section
			{
				for ( j = 0; j < size; j++)
					r[j] = r[j] - alpha * Akp[j] + alphaimg * Akpimg[j];
			}
			#pragma omp section
			{
				for ( j = 0; j < size; j++)
					rimg[j] = rimg[j] - alphaimg * Akp[j] - alpha * Akpimg[j];
			}
			#pragma omp section
			{
				for ( j = 0; j < size; j++)
					rtilde[j] = rtilde[j] - alpha * Akptilde[j] - alphaimg * Akptildeimg[j];
			}
			#pragma omp section
			{
				for ( j = 0; j < size; j++)
					rtildeimg[j] = rtildeimg[j] + alphaimg * Akptilde[j] - alpha * Akptildeimg[j];
			}
		}
			
//			if ( (p[j] != p[j])|| (p[j] + 1.0 ==p[j]) || (alpha != alpha)  || (alphaimg != alphaimg) || \
//				(pimg[j] != pimg[j]) )
//			{
//				printf("j = %d, r = (%lf, %lf), p = (%lf, %lf),alpha = (%lf, %lf).\n", \
//				j, r[j], rimg[j], p[j], pimg[j], alpha, alphaimg);
//				exit(0);
//			}


		//	BEGIN beta = (r^,r)_k+1/(r^, r)k		//
		comptempc = innerrold;
		comptempd = innerroldimg;	
		comptempa = omp_inner(size, r, rtilde, chunk) + omp_inner(size, rimg, rtildeimg, chunk);
		comptempb = omp_inner(size, rimg, rtilde, chunk) - omp_inner(size, r, rtildeimg, chunk);
		comptempnum = comptempa*comptempc + comptempb*comptempd;
		comptempnumimg = comptempb*comptempc - comptempa*comptempd;
		comptempden = comptempc*comptempc + comptempd*comptempd;
		beta = comptempnum/comptempden;
		betaimg = comptempnumimg/comptempden;

		innerrold = comptempa;
		innerroldimg = comptempb;

		//	END a = rho/ (r^, v)	//

		#pragma omp parallel for private(j) schedule(static, chunk)
		for ( j = 0; j < size; j++)		
		{
			pold[j] = p[j];
			poldimg[j] = pimg[j];
			ptildeold[j] = ptilde[j];
			ptildeoldimg[j] = ptildeimg[j];
		}

		#pragma omp parallel for private(j) schedule(static, chunk)
		for ( j = 0; j < size; j++)
		{
			p[j] = r[j] + beta * pold[j] - betaimg * poldimg[j];		//	p = r + beta * p	//
			pimg[j] = rimg[j] + betaimg * pold[j] + beta * poldimg[j];

			ptilde[j] = rtilde[j] + beta * ptildeold[j] + betaimg * ptildeoldimg[j];
			ptildeimg[j] = rtildeimg[j] - betaimg * ptildeold[j] + beta * ptildeoldimg[j];


			if ((p[j] !=p[j])||(p[j] + 1 ==p[j])||(pimg[j] !=pimg[j])||(pimg[j]+1==pimg[j]))
			{
				printf("s wrong, j = %d, alpha = %lf, r[j] = %lf, rimg[j] = %lf.\n", \
				j, alpha, r[j], rimg[j]);
				exit(0);
			}
		}



		preerrornorm = errornorm;
		errornorm = sqrt(omp_inner(size, r, r, chunk) + omp_inner(size, rimg, rimg, chunk))/normb;
//		printf("Run %d, error = %.8E\n", i, errornorm);


		loop_endtime = time(NULL);
		seconds = loop_endtime - loop_initime;
//		if (i == 1)
//			printf("This run took %ld minutes and %ld seconds.\n", seconds/60, seconds % 60);
		countrun++;

		if (errornorm < accuracy)
			break;
		if (errornorm != errornorm)
			break;
		if (fabs(preerrornorm - errornorm) < pow(10, -3) * accuracy)
			break;
		
	}

	printf("After %d runs, the error is %.10E.\n", countrun, errornorm);
	
	#pragma omp parallel sections private(i)
	{
		#pragma omp section
		{
			for ( i = 0; i < size; i++)
			{
				fprintf(fp3Dx, "%.15lf ", x[i]);
				if (x[i] != x[i])
				{
					printf("At i = %d, x is nan.\n", i);
					exit(0);
				}
			}
		}
		#pragma omp section
		{
			for ( i = 0; i < size; i++)
			{
				fprintf(fp3Dximg, "%.15lf ", ximg[i]);
				if (ximg[i] != ximg[i])
				{
					printf("At i = %d, x is nan.\n", i);
					exit(0);
				}
			}
		}
	}
	#pragma omp parallel for private(i) schedule(static, chunk)
	for(i = 0; i < size; i++)
	{
		Density[i] = x[i];
		Densityimg[i] = ximg[i];
	}

	fclose(fp3Dx);
	fclose(fp3Dximg);
//	fclose(fp3DC);
//	fclose(fp3DCimg);
	return 0;
}






int BiCGSTABcomplex(int size, double epsilon)
{
	int i, j, k, counter;
	FILE *fp3Dx, *fprealb, *fp3Dximg, *fpimgb;
//	FILE *fp3DC;		//	DEBUG	//	
//	FILE *fp3DCimg;		//	DEBUG	//

	int compsize = 2 * size;
	double nonhomogeneous = 0.0;
	double nonhomogeneousimg = 0.0;
	double zx, zy;
	double zstarx[size], zstary[size];
//	double epsilon = FACTOR * h;
	int count = 0;;

	double tempC[2];
	double r[size], rhat[size], b[size], x[size], v[size], p[size], s[size], t[size];
	double rimg[size], rhatimg[size], bimg[size], ximg[size], vimg[size], pimg[size], simg[size], timg[size];
	double error[size], errorimg[size];
	double alpha = 1.0, rhoi = 1.0, rhoi1 = 1.0, omega = 1.0;
	double alphaimg = 0.0, rhoiimg = 0.0, rhoi1img = 0.0, omegaimg = 0.0;
	double beta, betaimg;
	double errornorm, preerrornorm;
	double comptemp, comptempnum, comptempnumimg, comptempden, comptempa, comptempb, comptempc, comptempd;

	double sum = 0.0, regsum = 0.0, sumimg = 0.0, regsumimg = 0.0;
	double normb = 0.0;
	double BC_sum = 0.0;		//	for zero BC	//
	double accpower;
	double accuracy = DEFAULTACC;

	int countrun = 0;
	int checksize, checksizeimg;

	char denfilename[50], denimgfilename[50];

//	static int flag = 0;

	time_t loop_initime, loop_endtime, seconds;

	sprintf(denfilename, "DensityN%dk%d.txt", (int)(N), (int)(WAVE));
	if ((fp3Dx = fopen(denfilename, "w+"))==NULL)
	{
		printf("Open file %s error.\n", denfilename);
		return 0;
	}
	sprintf(denimgfilename, "DensityimgN%dk%d.txt", (int)(N), (int)(WAVE));
	if ((fp3Dximg = fopen(denimgfilename, "w+"))==NULL)
	{
		printf("Open file %s error.\n", denimgfilename);
		return 0;
	}
//	fprintf(fp3Dx, "%d\n", size);

	if ((fprealb = fopen("largerealb.txt", "r"))==NULL)
	{
		printf("Open file large3Drealb.txt error.\n");
		return 0;
	}
	if ((fpimgb = fopen("largeimgb.txt", "r"))==NULL)
	{
		printf("Open file largeimgb.txt error.\n");
		return 0;
	}
//	if ((fp3DC = fopen("large3DC.txt", "w+"))==NULL)
//	{
//		printf("Open file large3DC.txt error.\n");
//		return 0;
//	}
//	fprintf(fp3DC, "%d\n", size);
//	if ((fp3DCimg = fopen("large3DCimg.txt", "w+"))==NULL)
//	{
//		printf("Open file large3DCimg.txt error.\n");
//		return 0;
//	}
//	fprintf(fp3DCimg, "%d\n", size);

	printf("Enter the accuracy in 10^x ( 0 > x > -10)\n");
	scanf("%lf", &accpower);
	if ( (accpower > -10.0) && (accpower < 0.0) )
		accuracy = pow(10, accpower);
	printf("Using accuracy = %.8E.\n", accuracy);

	//	INITIALIZING THE THINGS NEEDED FOR BICONJUGATE GRADIENT STABILIZED	//

	fscanf(fprealb, "%d", &checksize);
	fscanf(fpimgb, "%d", &checksizeimg);
	if ((checksize != size) || (checksizeimg!=size))
	{
		printf("size error. Realsize = %d, realpart = %d, imgpart = %d\n", size, checksize, checksizeimg);
		return 0;
	}
	for ( counter = 0; counter < size; counter++)
	{
		fscanf(fprealb, "%lf", &b[counter]);
		fscanf(fpimgb, "%lf", &bimg[counter]);

		if ((b[counter]!=b[counter])||(b[size+counter]!=b[size+counter]))
		{
			printf("b[%d] = %lf.\n", counter, b[counter]);
			printf("bimg[%d] = %lf.\n", counter, bimg[counter]);
			exit(0);
		}
	}
	normb = sqrt(inner(size, b, b) + inner(size, bimg, bimg));
	for ( counter = 0; counter < size; counter++)
	{
		r[counter] = b[counter];
		rimg[counter] = bimg[counter];
		v[counter] = 0.0;
		vimg[counter] = 0.0;
		p[counter] = 0.0;
		pimg[counter] = 0.0;
		rhat[counter] = r[counter];
		rhatimg[counter] = rimg[counter];
		if ( (rhat[counter] != r[counter]) || (rhatimg[counter] != rimg[counter]) )
		{
			printf("counter = %d.\n, rhat = %lf, r = %lf.\n", counter, rhat[counter], r[counter]);
			printf("rhatimg = %lf, rimg = %lf.\n", rhatimg[counter], rimg[counter]);
		}
	}
	fclose(fprealb);
	fclose(fpimgb);
	///////////////////////////////////////////////////////////////////////////////////////


	for ( i = 0; i < size; i++)
	{
		BC_sum += fabs(r[i]);
		BC_sum += fabs(rimg[i]);
	}
	if ( BC_sum < pow(10,-13.0))	//	for zero BC	//
	{
		printf("Here.\n");
		r[0] = 0.000000001;
		rimg[0] = 0.000000001;
		rhat[0] = 0.000000001;
		rhatimg[0] = 0.000000001;
	}

	if ( BC_sum < pow(10, -15.0) )
	{
		printf("Zero BC.\n");
		for ( i = 0; i < size; i++)
		{
			fprintf(fp3Dx, "%.15lf ", 0.0);
			fprintf(fp3Dximg, "%.15lf ", 0.0);
		}
		fclose(fp3Dx);
		fclose(fp3Dximg);
		return 0;
	}


	for ( i = 1; i <= 20000; i++)
	{
		loop_initime = time(NULL);

		rhoi1 = rhoi;					//	rhoi-1 = rhoi	//
		rhoi1img = rhoiimg;

//		rhoi = inner(size, rhat, r) + inner(size, rhatimg, rimg);			//	rho = (r^, r)			//		
//		rhoiimg = inner(size, rhatimg, r) - inner(size, rhat, rimg);
		rhoi = inner(size, rhat, r) - inner(size, rhatimg, rimg);
		rhoiimg = inner(size, rhatimg, r) + inner(size, rhat, rimg);

		//	DEFINE IF A = a+bi, B = c+di, then (A,B) = sigma(A*B') = (a+bi)dot(c-di)	//
		//	THAT IS, CONJUGATE THE LATTER VECTOR	//

//		printf("rhoi = %lf, rhat[0] = %lf, r[0] = %lf.\n", rhoi, rhat[0], r[0]);
//		if (rhoi!=rhoi)
//			printf("rhoi.\n");

		//	BEGIN beta = (rhoi*alpha)/(rhoi1*omega)	//
		comptempa = rhoi*alpha - rhoiimg*alphaimg;			
		comptempb = rhoiimg*alpha + rhoi*alphaimg;
		comptempc = rhoi1*omega - rhoi1img*omegaimg;
		comptempd = rhoi1*omegaimg + rhoi1img*omega; 
		comptempnum = (comptempa*comptempc + comptempb*comptempd);
		comptempnumimg = ( comptempb*comptempc - comptempa*comptempd );
		comptempden = (comptempc*comptempc + comptempd*comptempd);
		beta = comptempnum/comptempden;				//	beta = (rhoi/rhoi1) * (a/w)	//
		betaimg = comptempnumimg/comptempden;
		//		printf("beta = %lf\n", beta);
		//	END beta = (rhoi*alpha)/(rhoi1*omega)	//


		for ( j = 0; j < size; j++)
		{
			p[j] = r[j] + beta * (p[j] - omega*v[j] + omegaimg*vimg[j]) - betaimg*(pimg[j] - omegaimg*v[j] - omega*vimg[j]);
			pimg[j] = rimg[j] + betaimg*(p[j] - omega*v[j] + omegaimg*vimg[j]) + beta*(pimg[j] - omegaimg*v[j] - omega*vimg[j]);
			//	p = r + beta( p - w * v )  //

//			if ( (p[j] != p[j])|| (p[j] + 1.0 ==p[j]) )
//				printf("j = %d, r[j] = %lf, v[j] = %lf, beta = %lf.\n", j, r[j], v[j], beta);
		}

		//////	FIRST MATRIX MULTIPLICATION	////////////
		for ( j = 0; j < size; j++)
		{
			sum = 0.0;
			sumimg = 0.0;
			for ( k = 0; k < size; k++)
			{
//				if (i==1)
//				{
//					fprintf(fp3DC, "%.15lf ", get_C(size, j, k));
//				}
				get_Ccomp(size, j, k, epsilon, tempC);
				sum += (tempC[0] * p[k] - tempC[1] * pimg[k]);
				sumimg += (tempC[0] * pimg[k] + tempC[1] * p[k]);
			}
//			if (i==1)
//				fprintf(fp3DC, "\n");
			v[j] = sum;				//	v = A * p	//
			vimg[j] = sumimg;

			if ((v[j] != v[j])||(v[j] +1 == v[j])||(vimg[j] != vimg[j])||(vimg[j]+1==vimg[j]))
			{
				printf("v wrong, j = %d, v[j] = %lf, vimg[j] = %lf.\n", j, v[j], vimg[j]);
				exit(0);
			}
//			printf("v[%d] = %.15lf\n", j, v[j]);
		}
		loop_endtime = time(NULL);
		seconds = loop_endtime - loop_initime;
		printf("Made to 1st. Took %ld minutes %ld seconds.\n", seconds/60, seconds % 60);
		/////	FIRST MATRIX MULTIPLICATION ENDS	////

//		printf("v[0] = %.15lf.\n", v[0]);
//		printf("rhoi = %.15lf.\n", rhoi);

		//	BEGIN a = rho/(r^, v)		//
//		comptempc = inner(size, rhat, v) + inner(size, rhatimg, vimg);
//		comptempd = inner(size, rhatimg, v) - inner(size, rhat, vimg);
		comptempc = inner(size, rhat, v) - inner(size, rhatimg, vimg);
		comptempd = inner(size, rhatimg, v) + inner(size, rhat, vimg);
		comptempnum = rhoi*comptempc + rhoiimg*comptempd;
		comptempnumimg = rhoiimg*comptempc - rhoi*comptempd;
		comptempden = comptempc*comptempc + comptempd*comptempd;
		alpha = comptempnum/comptempden;
		alphaimg = comptempnumimg/comptempden;
		//	END a = rho/ (r^, v)	//

		for ( j = 0; j < size; j++)		
		{
			s[j] = r[j] - alpha * v[j] + alphaimg * vimg[j];	//	s = r - a * v	//
			simg[j] = rimg[j] - alphaimg * v[j] - alpha * vimg[j];
			if ((s[j] !=s[j])||(s[j] + 1 ==s[j])||(simg[j] !=simg[j])||(simg[j]+1==simg[j]))
			{
				printf("s wrong, j = %d, alpha = %lf, r[j] = %lf, rimg[j] = %lf.\n", j, alpha, r[j], rimg[j]);
				exit(0);
			}
		}

		//////	SECOND MATRIX MULTIPLICATION	/////////////
		for ( j = 0; j < size; j++)
		{
			sum = 0.0;
			sumimg = 0.0;
			for ( k = 0; k < size; k++)
			{
				get_Ccomp(size, j, k, epsilon, tempC);
				sum += (tempC[0]*s[k] - tempC[1]*simg[k]);
				sumimg += (tempC[0]*simg[k] + tempC[1]*s[k]);
			}
			t[j] = sum;				//	t = A * s		//
			timg[j] = sumimg;
			if ((t[j]!=t[j])||(timg[j]!=timg[j]))
			{
				printf("t wrong, j = %d, t[j] = %lf, timg[j] = %lf.\n", j, t[j], timg[j]);
				printf("alpha = %lf, r[%d] = %lf.\n", alpha, j, r[j]);
				exit(0);
			}
		}
		///////	SECOND MATRIX MULTIPLICATION ENDS	//////

		//	BEGIN	w = (t, s)/(t, t)	//
//		comptempa = inner(size, t, s) + inner(size, timg, simg);
//		comptempb = inner(size, timg, s) - inner(size, t, simg);
		comptempa = inner(size, t, s) - inner(size, timg, simg);
		comptempb = inner(size, timg, s) + inner(size, t, simg);
		comptempden = inner(size, t, t) + inner(size, timg, timg);
		omega = comptempa/comptempden;
		omegaimg = comptempb/comptempden;
		//	END 	w = (t, s)/(t, t)	//

		for ( j = 0; j < size; j++)
		{
			x[j] = x[j] + alpha*p[j] - alphaimg*pimg[j] + omega*s[j] - omegaimg*simg[j];	//	x = x + a * p + w * s	//
			ximg[j] = ximg[j] + alphaimg*p[j] + alpha*pimg[j] + omegaimg*s[j] + omega*simg[j];
		}
	
		///////	THIRD MATRIX MULTIPLICATION	/////////////
//		for (int j = 0; j < size; j++)
//		{
//			sum = 0.0;
//			for (int k = 0; k < size; k++)
//			{
//				sum += (get_C(size, j, k) * x[k]);
//			}
//			error[j] = sum - b[j];			//	err = b - A * x		//
//			if (error[j] != error[j])
//				printf("error wrong, j = %d, x = %lf, b = %lf, omega = %lf.", j, x[j], b[j], omega);
//		}

		for ( j = 0; j < size; j++)
		{
			r[j] = s[j] - omega * t[j] + omegaimg * timg[j];		//	r = s - w * t		//
			rimg[j] = simg[j] - omegaimg * t[j] - omega * timg[j];
		}

		preerrornorm = errornorm;
		errornorm = sqrt(inner(size, r, r) + inner(size, rimg, rimg))/normb;
		printf("Run %d, error = %.8E\n", i, errornorm);


		loop_endtime = time(NULL);
		seconds = loop_endtime - loop_initime;
		printf("This run took %ld minutes and %ld seconds.\n", seconds/60, seconds % 60);
		countrun++;

		if (errornorm < accuracy)
			break;
		if (errornorm != errornorm)
			break;
		if (fabs(preerrornorm - errornorm) < pow(10, -3) * accuracy)
			break;
		
	}

	printf("After %d runs, the error is %.10E.\n", countrun, errornorm);
	
	for ( i = 0; i < size; i++)
	{
		fprintf(fp3Dx, "%.15lf ", x[i]);
		fprintf(fp3Dximg, "%.15lf ", ximg[i]);
	}
	
	fclose(fp3Dx);
	fclose(fp3Dximg);
//	fclose(fp3DC);
	return 0;
}


double get_C(int size, int pi, int pj)
{
	FILE *fp, *fp2;

	int i, j;
	int counti = 0;
	int countj = 0;

	static double *zx, *zy, *wx, *wy;
	static double *zstarx, *zstary, *wstarx, *wstary;
	static double *J, *delta;
	static double *gradx, *grady;
	static double *result;
	static double epsilon;
	static double tau, threshold, regfactor, regfactorimg;
	
	static double **matrix;


	double tempzx, tempzy, tstarx, tstary;
	double starxi, staryi, starxj, staryj;
	double C, C_reg;
	double partial, partial_reg;
	double partialx, partialy;
	double thisdist, thisgradx, thisgrady;
//	double h = (INT_U - INT_L)/N;
//	double epsilon = FACTOR * h;
//	double tau = h;
//	time_t begin, middle, end, seconds;

	static int ini = 0;
	static int flag = 0;
	static unsigned int fastflag = 0;
	static unsigned int fillinflag = 0;
	static long int debugcount = 0;

	if (size <= FASTSIZEBOUND)
	{
		if ( fastflag >= size*size )
			return matrix[pi][pj];
		fastflag++;
	}

	//	TRY TO SAVE MORE TIME	//
//	else if ( ( pi < FASTSIZEBOUND) && (pj < FASTSIZEBOUND) && (fillinflag >= FASTSIZEBOUND*FASTSIZEBOUND) )
//	{
//		return matrix[pi][pj];
//	}


	if (ini == 0)	//	FIRST TIME, THE ONLY TIME DOING INITIALIZATION FOR X, Y, Z POINTS //
	{
		tau = H/TAU_FACTOR;
		threshold = tau * tau;
		epsilon = FACTOR * H;

		regfactor = 1.0/(4.0*PI*radius);

		zx = (double *) malloc(size * sizeof(double));
		zy = (double *) malloc(size * sizeof(double));
		zstarx = (double *) malloc(size * sizeof(double));
		zstary = (double *) malloc(size * sizeof(double));
		delta = (double *) malloc(size * sizeof(double));
		J = (double *) malloc(size * sizeof(double));
		gradx = (double *) malloc(size * sizeof(double));
		grady = (double *) malloc(size * sizeof(double));
		result= (double *) malloc(size * sizeof(double));
		

		if ( (zx==NULL)||(zy==NULL) || (zstarx==NULL) || (zstary==NULL) )
			printf("Memory allocation error....points.\n");
		if ( (delta==NULL)||(J==NULL)||(gradx==NULL)|| (grady==NULL) )
			printf("Memory allocation error....calcs.\n");
		if (result==NULL)
			printf("Memory allocation error....result.\n");

		if (ABS_MODE != 1)
			epsilon = FACTOR * H;

		for ( i = 0; i < SIZE; i++)
		{
			tempzy = INT_L + H * (double)(i);
			for ( j = 0; j < SIZE; j++)
			{
				tempzx = INT_L + H * (double)(j);
				if (ABS_MODE == 1)
				{
					epsilon = ABS_EPSILON;
				}
				else if (ABS_MODE == 0)
				{
					thisgradx = gradientdx(tempzx, tempzy);
					thisgrady = gradientdy(tempzx, tempzy);
//					epsilon = ( fabs(thisgradx) + fabs(thisgrady) ) * FACTOR * H;
				}
				//	if not in tubular region, don't calculate.	//
//				thisdist = distance(tempzx, tempzy);
				thisdist = DIST[i][j];
				if (fabs(thisdist) > epsilon)
					continue;

				zx[counti] = tempzx;
				zy[counti] = tempzy;
//				tstarx = starx(tempzx, tempzy);
//				tstary = stary(tempzx, tempzy);
				tstarx = tempzx - thisdist * DIST_DX[i][j];
				tstary = tempzy - thisdist * DIST_DY[i][j];

				zstarx[counti] = tstarx;
				zstary[counti] = tstary;

				delta[counti] = cal_delta(epsilon, thisdist);
//				J[counti] = cal_J(J_M, tempzx, tempzy);
				J[counti] = 1. + thisdist * CURV[i][j];
//				gradx[counti] = cal_interpgrad(1, tstarx, tstary);
//				grady[counti] = cal_interpgrad(2, tstarx, tstary);
				cal_interpgrad(tstarx, tstary, DIST_DX, DIST_DY, &gradx[counti], &grady[counti]);

				result[counti] = H*H * delta[counti] * J[counti];

				counti++;	// in the end should total up to size	//
			}
		}
		
		if (size <= FASTSIZEBOUND)
		{
			matrix = (double **) malloc(size * sizeof(double *));
			for ( i = 0; i < size; i++)
			{
				matrix[i] = (double *) malloc(size * sizeof(double));
				if (matrix[i] == NULL)
				{
					printf("Memory allocation error for matrix[%d].\n", i);
					exit(0);
				}
			}
		}
		else
		{
			matrix = (double **) malloc(FASTSIZEBOUND * sizeof(double *));
			for ( i = 0; i < FASTSIZEBOUND; i++)
			{
				matrix[i] = (double *) malloc(FASTSIZEBOUND *sizeof(double));
				if (matrix[i] == NULL)
				{
					printf("Memory allocation error for matrix[%d].\n", i);
					exit(0);
				}
			}
		}

		ini = 1;	//	DENOTING DATA INITIALIZED	//
	}


	starxi = zstarx[pi];
	staryi = zstary[pi];
	starxj = zstarx[pj];
	staryj = zstary[pj];

	partialx = cal_partial(1, starxi, staryi, starxj, staryj);// x coordinate	//
	partialy = cal_partial(2, starxi, staryi, starxj, staryj);// y coordinate	//

	//	for interior problem	//
	if (INTERIOR_MODE==1)
		partial = (-1.0) * ( partialx * gradx[pj] + partialy * grady[pj] );	//	outer normal
	//	for exterior problem	//
	else
		partial = partialx * gradx[pj] + partialy * grady[pj];			//	inner normal


	if ( ( (starxi-starxj)*(starxi-starxj) + (staryi-staryj)*(staryi-staryj) )< threshold )
	{
		partial_reg = regfactor;
	}
	else
	{
		partial_reg = partial;
	}

	//	calculate matrix elements	//

	C = result[pj] * partial ;		//	( h^3*deltaj*Jj ) * partial	//
	C_reg = result[pj] * partial_reg;

	//	C = B + 1/2 I	//
	if (pi == pj)
	{
		C += 0.5;	//(delta*cal_J(J_M, wx, wy));
		C_reg += 0.5;	//(delta*cal_J(J_M, wx, wy));
//		C += 1.0;	//	+1 is for Helmholtz, + 0.5 is for Laplace->which is wrong.
//		C_reg += 1.0;	//	+1 is for Helmholtz, + 0.5 is for Laplace. Both should be + 0.5
	}

	if (size <= FASTSIZEBOUND)
	{
		if (REG_MODE==0)
			matrix[pi][pj] = C;
		else
			matrix[pi][pj] = C_reg;
	}
	//	TRY TO SAVE MORE	//
//	else if (fillinflag < FASTSIZEBOUND*FASTSIZEBOUND)
//	{
//		if ( (pi < FASTSIZEBOUND) && (pj < FASTSIZEBOUND) )
//		{
//			if (REG_MODE==0)
//				matrix[pi][pj] = C;
//			else
//				matrix[pi][pj] = C_reg;
//			fillinflag++;
//		}
//	}

	if ((C_reg!=C_reg)||(C_reg+1==C_reg))
	{
		printf("i = %d, j = %d, partial = %.10lf, C_reg = %.10lf.\n", pi,  pj, partial, C_reg);
		exit(0);
	}
	if ((C!=C)||(C+1==C))
	{
		printf("i = %d, j = %d, partial = %.10lf, C = %.10lf.\n", pi, pj, partial_reg, C);
		exit(0);
	}

	if (REG_MODE==0)
		return C;
	else
		return C_reg;
}


double get_Ccomp(int size, int pi, int pj, double epsilon, double tempC[])
{
	int i, j;
	static FILE *fp, *fp2;
//	static FILE *fpxstar, *fpystar;

	int counti = 0;
	int countj = 0;

	static double *zx, *zy, *wx, *wy;
	static double *zstarx, *zstary, *wstarx, *wstary;
	static double *J, *delta;
	static double *gradx, *grady;
	static double *result, *regfactor;
	static double tau, threshold;

	static double **matrix;
	static double **matriximg;
	static double *temp;


	double tempzx, tempzy, tstarx, tstary;
	double starxi, staryi, starxj, staryj;
	double C, C_reg;
	double partial, partial_reg, partialimg, partialimg_reg;
	double partialx, partialy, partialximg, partialyimg;
	double thisdist, thisgradx, thisgrady, thiscurv;

	double parameter, commonterm, thisnorm, xminusy1, xminusy2, xmyinner_pd, gradnorm;
	double BesselY1value, BesselJ1value;

//	double H = (INT_U - INT_L)/N;
//	double epsilon = FACTOR * H;
//	double tau = H;
//	time_t begin, middle, end, seconds;

	static int ini = 0;
	static int flag = 0;
	static unsigned long int fastflag = 0;
	static unsigned long int fillinflag = 0;


	if ( (pi == -1) && (pj == -1) )
	{
		for ( i = 0; i < size; i++)
		{
			free(matrix[i]);
			free(matriximg[i]);
		}
		free(matrix);
		free(matriximg);
		free(zx);
		free(zy);
		free(zstarx);
		free(zstary);
		free(delta);
		free(J);
		free(gradx);
		free(grady);
		free(result);
		free(regfactor);
		return 0.;
	}

	if (size <= FASTSIZEBOUND)
	{
		if (fastflag >= size*size)
		{
			tempC[0] = matrix[pi][pj];
			tempC[1] = matriximg[pi][pj];
			return 0;
		}
		fastflag++;
	}
	else if ((pi < FASTSIZEBOUND) && (pj < FASTSIZEBOUND) && (fillinflag >= (unsigned long int)FASTSIZEBOUND * (unsigned long int )FASTSIZEBOUND))
	{
		tempC[0] = matrix[pi][pj];
		tempC[1] = matriximg[pi][pj];
		return 0;
	}

	if (ini == 0)	//	FIRST TIME, THE ONLY TIME DOING INITIALIZATION FOR X, Y, Z POINTS //
	{
		tau = H/TAU_FACTOR;
		threshold = tau * tau;
//		epsilon = FACTOR * H;

//		if (INTERIOR_MODE==1)
//			regfactor = 1.0/(4.0*PI*radius);
//		else
//			regfactor = -1.0/(4.0*PI*radius);
//		regfactorimg = 0.0;

		zx = (double *) malloc(size * sizeof(double));
		zy = (double *) malloc(size * sizeof(double));
		zstarx = (double *) malloc(size * sizeof(double));
		zstary = (double *) malloc(size * sizeof(double));
		delta = (double *) malloc(size * sizeof(double));
		J = (double *) malloc(size * sizeof(double));
		gradx = (double *) malloc(size * sizeof(double));
		grady = (double *) malloc(size * sizeof(double));
		result= (double *) malloc(size * sizeof(double));
		regfactor = (double *) malloc(size * sizeof(double));
		

		if ( (zx==NULL)||(zy==NULL) || (zstarx==NULL) || (zstary==NULL) )
			printf("Memory allocation error....points.\n");
		if ( (delta==NULL)||(J==NULL)||(gradx==NULL)|| (grady==NULL) )
			printf("Memory allocation error....calcs.\n");
		if ((result==NULL) || (regfactor==NULL) )
			printf("Memory allocation error....result.\n");

//		if ( (fpxstar = fopen("xstar.txt", "w+")) == NULL )
//		{
//			printf("Open file xstar.txt error.\n");
//			exit(0);
//		}
//		if ( (fpystar = fopen("ystar.txt", "w+")) == NULL )
//		{
///			printf("open file ystar.txt error.\n");
//			exit(0);
//		}

//		if (ABS_MODE != 1)
//			epsilon = FACTOR * H;


		for ( i = 0; i < SIZE; i++)
		{
			tempzy = INT_L + H * (double)(i);
			for ( j = 0; j < SIZE; j++)
			{
				tempzx = INT_L + H * (double)(j);

//				if (ABS_MODE == 1)
//				{
//					epsilon = ABS_EPSILON;
//				}
//				else if (ABS_MODE == 0)
//				{
//					thisgradx = gradientdx(tempzx, tempzy);
//					thisgrady = gradientdy(tempzx, tempzy);
//					epsilon = ( fabs(thisgradx) + fabs(thisgrady) ) * FACTOR * h;
//				}
//				thisdist = distance(tempzx, tempzy);
				thisdist = DIST[i][j];
				//	if not in tubular region, don't calculate.	//

				if (TRAD_ONESIDED == 0)
				{
					if ( fabs(thisdist) > epsilon ) 
						continue;
				}
				else
				{
					if ( (thisdist > epsilon ) || (thisdist < 0.) )
						continue;
				}

				zx[counti] = tempzx;
				zy[counti] = tempzy;
//				tstarx = starx(tempzx, tempzy);
//				tstary = stary(tempzx, tempzy);
				tstarx = tempzx - thisdist * DIST_DX[i][j];
				tstary = tempzy - thisdist * DIST_DY[i][j];
				

				zstarx[counti] = tstarx;
				zstary[counti] = tstary;
//				fprintf(fpxstar, "%lf ", tstarx);
//				fprintf(fpystar, "%lf ", tstary);

				if (TRAD_ONESIDED == 0)
					delta[counti] = cal_delta(epsilon, thisdist);
				else
					delta[counti] = 2. * cal_delta(epsilon, thisdist);

//				J[counti] = cal_J(J_M, tempzx, tempzy);
				J[counti] = 1. + thisdist * CURV[i][j];
				thiscurv = cal_curv(tstarx, tstary, DIST);
				regfactor[counti] = thiscurv/(4.*PI);
//				gradx[counti] = cal_interpgrad(1, tstarx, tstary);
//				grady[counti] = cal_interpgrad(2, tstarx, tstary);
				cal_interpgrad(tstarx, tstary, DIST_DX, DIST_DY, &gradx[counti], &grady[counti]);
				gradnorm = sqrt(gradx[counti] * gradx[counti] + grady[counti]*grady[counti]);
				if (fabs(gradnorm - 1.0) > 0.1)
				{
					printf("At (%lf, %lf), grad (%lf, %lf), norm %lf.\n", \
					tstarx, tstary, gradx[counti], grady[counti], gradnorm);
					exit(0);
				}

				result[counti] = HSQUARED * delta[counti] * J[counti];

				counti++;	// in the end should total up to size	//
			}
		}

//		fclose(fpxstar);
//		fclose(fpystar);

		if (size <= FASTSIZEBOUND)
		{
			matrix = (double **) malloc(size * sizeof(double *));
			matriximg = (double **) malloc(size * sizeof(double *));
			for ( i = 0; i < size; i++)
			{
				matrix[i] = (double *) malloc(size * sizeof(double));
				matriximg[i] = (double *) malloc(size * sizeof(double));
				if ( (matrix[i]==NULL) || (matriximg[i]==NULL) )
				{
					printf("Memory allocation error for matrix[%d].\n", i);
					exit(0);
				}
			}
		}
		else
		{
			matrix = (double **) malloc(FASTSIZEBOUND * sizeof(double *));
			matriximg = (double **) malloc(FASTSIZEBOUND * sizeof(double *));
			for ( i = 0; i < FASTSIZEBOUND; i++)
			{
				matrix[i] = (double *) malloc(FASTSIZEBOUND * sizeof(double));
				matriximg[i] = (double *) malloc(FASTSIZEBOUND * sizeof(double));
				if ( (matrix[i]==NULL) || (matriximg[i]==NULL) )
				{
					printf("Memory allocation error for matrix[%d].\n", i);
					exit(0);
				}
			}
		}

		ini = 1;	//	DENOTING DATA INITIALIZED	//
	}


	starxi = zstarx[pi];
	staryi = zstary[pi];
	starxj = zstarx[pj];
	staryj = zstary[pj];


//	if (DIRICHLET == 1)
//	{
///////////////////////////////////////////////		faster kernel		//////////////////////////////

	xminusy1 = starxi-starxj;
	xminusy2 = staryi-staryj;

	thisnorm = sqrt( xminusy1*xminusy1 + xminusy2*xminusy2 );


	if (thisnorm > tau )
	{

		parameter= WAVE * thisnorm;
		BesselY1value = BesselY1(parameter);
		BesselJ1value = BesselJ1(parameter);

		commonterm = 0.25 * WAVE / thisnorm;

//		partialx = commonterm * BesselY1value * xminusy1;
//		partialy =  commonterm * BesselY1value * xminusy2;
//		partialximg = -1. * commonterm * BesselJ1value * xminusy1;
//		partialyimg = -1. * commonterm * BesselJ1value * xminusy2;	

		if (DIRICHLET == 1)
			xmyinner_pd = xminusy1 * gradx[pj] + xminusy2 * grady[pj];
		else
			xmyinner_pd = xminusy1 * gradx[pi] + xminusy2 * grady[pi];

/////////////////////////////////////////////////////////////////////////////////////////////////////////////


//		partialx = cal_partial(1, starxi, staryi, starxj, staryj);// x coordinate real	//
//		partialy = cal_partial(2, starxi, staryi, starxj, staryj);// y coordinate real	//
//		partialximg = cal_partial(3, starxi, staryi, starxj, staryj);// x coordinate imaginary //
//		partialyimg = cal_partial(4, starxi, staryi, starxj, staryj);// y coordinate imaginary //
//	}
//	else
//	{
//		partialx = cal_partial(1, starxj, staryj, starxi, staryi);// x coordinate real	//
//		partialy = cal_partial(2, starxj, staryj, starxi, staryi);// y coordinate real	//
//		partialximg = cal_partial(3, starxj, staryj, starxi, staryi);// x coordinate imaginary //
//		partialyimg = cal_partial(4, starxj, staryj, starxi, staryi);// y coordinate imaginary //

//		partialx = cal_partial(1, starxi, staryi, starxj, staryj);// x coordinate real	//
//		partialy = cal_partial(2, starxi, staryi, starxj, staryj);// y coordinate real	//
//		partialximg = cal_partial(3, starxi, staryi, starxj, staryj);// x coordinate imaginary //
//		partialyimg = cal_partial(4, starxi, staryi, starxj, staryj);// y coordinate imaginary //
//	}


		//	for interior problem	//
		if (INTERIOR_MODE==1)
		{
			partial = (-1.0) * commonterm * BesselY1value * xmyinner_pd;	//	outer normal
			partialimg = commonterm * BesselJ1value * xmyinner_pd;
		}
		//	for exterior problem	//
		else
		{
			partial =-1.0 *  commonterm * BesselY1value * xmyinner_pd;	//	outer normal
			partialimg = commonterm * BesselJ1value * xmyinner_pd;
		}
	}
//	if ( ( (starxi-starxj)*(starxi-starxj) + (staryi-staryj)*(staryi-staryj) )< threshold )
//	{
//		partial_reg = regfactor;
//		partialimg_reg = regfactorimg;
//	}
	else
	{
		partial = regfactor[pi];
		partialimg_reg = 0.0;
//		partial_reg = partial;
//		partialimg_reg = partialimg;
	}

	//	calculate matrix elements	//

	tempC[0] = result[pj] * partial ;		//	( h^3*deltaj*Jj ) * partial	//
	tempC[1] = result[pj] * partialimg;
//	if (REG_MODE==1)
//	{
//		tempC[0] = result[pj] * partial_reg;
//		tempC[1] = result[pj] * partialimg_reg;
//	}

	//	C = B + 1/2 I	//
	if (pi == pj)
	{
		if (DIRICHLET == 1)
			tempC[0] += 0.5;	//(delta*cal_J(J_M, wx, wy));
		else
			tempC[0] -= 0.5;
//		if (REG_MODE==1)
//		{
//			C_reg += 0.5;	//(delta*cal_J(J_M, wx, wy));
//		}
//		C += 1.0;	//	+1 is for Helmholtz, + 0.5 is for Laplace->which is wrong.
//		C_reg += 1.0;	//	+1 is for Helmholtz, + 0.5 is for Laplace. Both should be + 0.5
	}

	if (size <= FASTSIZEBOUND)
	{
		matrix[pi][pj] = tempC[0];
		matriximg[pi][pj] = tempC[1];
	}
	else if ( (pi < FASTSIZEBOUND) && (pj < FASTSIZEBOUND) && (fillinflag < (unsigned long int)FASTSIZEBOUND*(unsigned long int)FASTSIZEBOUND) )
	{
		matrix[pi][pj] = tempC[0];
		matriximg[pi][pj] = tempC[1];
		fillinflag++;
	}

	return 0.;
//	if ((C_reg!=C_reg)||(C_reg+1==C_reg))
//	{
//		printf("i = %d, j = %d, partial = %.10lf, C_reg = %.10lf.\n", pi,  pj, partial, C_reg);
//		exit(0);
//	}
	if ((tempC[0]!=tempC[0])||(tempC[0]+1==tempC[0])||(tempC[1]!=tempC[1])||(tempC[1]+1==tempC[1]))
	{
		printf("i = %d, j = %d, partial = %lf, C = %.10lf, Cimg = %.10lf.\n", \
			pi, pj, partial, tempC[0], tempC[1]);
		exit(0);
	}



//	if (REG_MODE==0)
//		return C;
//	else
//		return C_reg;
}


double omp_get_Ccomp(int size, int pi, int pj, double epsilon, double tempC[])
{
	int i, j;
	static FILE *fp, *fp2;
//	static FILE *fpxstar, *fpystar;

	int counti = 0;
	int countj = 0;

	static double *zx, *zy, *wx, *wy;
	static double *zstarx, *zstary, *wstarx, *wstary;
	static double *J, *delta;
	static double *gradx, *grady, *regfactor;
	static double *result;
	static double tau, threshold;

	static double *temp;


	double tempzx, tempzy, tstarx, tstary;
	double starxi, staryi, starxj, staryj;
	double C, C_reg;
	double partial, partial_reg, partialimg, partialimg_reg;
	double partialx, partialy, partialximg, partialyimg;
	double thisdist, thisgradx, thisgrady, thiscurv;

	static int ini = 0;
	static int flag = 0;
//	static unsigned int fastflag = 0;
//	static unsigned int fillinflag = 0;
	double parameter, commonterm, thisnorm, xminusy1, xminusy2, xmyinner_pd, gradnorm;
	double BesselY1value, BesselJ1value;


	if ( (pi == -1) && (pj == -1) )
	{
		free(zx);
		free(zy);
		free(zstarx);
		free(zstary);
		free(delta);
		free(J);
		free(gradx);
		free(grady);
		free(result);
		free(regfactor);
		return 0.;
	}

	if (ini == 0)	//	FIRST TIME, THE ONLY TIME DOING INITIALIZATION FOR X, Y, Z POINTS //
	{
		tau = H/TAU_FACTOR;
		threshold = tau * tau;

//		if (INTERIOR_MODE==1)
//			regfactor = 1.0/(4.0*PI*radius);
//		else
//			regfactor = -1.0/(4.0*PI*radius);
//		regfactorimg = 0.0;

		zx = (double *) malloc(size * sizeof(double));
		zy = (double *) malloc(size * sizeof(double));
		zstarx = (double *) malloc(size * sizeof(double));
		zstary = (double *) malloc(size * sizeof(double));
		delta = (double *) malloc(size * sizeof(double));
		J = (double *) malloc(size * sizeof(double));
		gradx = (double *) malloc(size * sizeof(double));
		grady = (double *) malloc(size * sizeof(double));
		result= (double *) malloc(size * sizeof(double));
		regfactor = (double *) malloc(size * sizeof(double));
		

		if ( (zx==NULL)||(zy==NULL) || (zstarx==NULL) || (zstary==NULL) )
			printf("Memory allocation error....points.\n");
		if ( (delta==NULL)||(J==NULL)||(gradx==NULL)|| (grady==NULL) )
			printf("Memory allocation error....calcs.\n");
		if ( (result==NULL) || (regfactor == NULL) )
			printf("Memory allocation error....result.\n");


		for ( i = 0; i < SIZE; i++)
		{
			tempzy = INT_L + H * (double)(i);
			for ( j = 0; j < SIZE; j++)
			{
				tempzx = INT_L + H * (double)(j);

//				thisdist = distance(tempzx, tempzy);
				thisdist = DIST[i][j];
				//	if not in tubular region, don't calculate.	//

				if (TRAD_ONESIDED == 0)
				{
					if ( fabs(thisdist) > epsilon ) 
						continue;
				}
				else
				{
					if ( (thisdist > epsilon ) || (thisdist < 0.) )
						continue;
				}

				zx[counti] = tempzx;
				zy[counti] = tempzy;
//				tstarx = starx(tempzx, tempzy);
//				tstary = stary(tempzx, tempzy);

				tstarx = tempzx - thisdist * DIST_DX[i][j];
				tstary = tempzy - thisdist * DIST_DY[i][j];

				zstarx[counti] = tstarx;
				zstary[counti] = tstary;

				if (TRAD_ONESIDED == 0)
					delta[counti] = cal_delta(epsilon, thisdist);
				else
					delta[counti] = 2. * cal_delta(epsilon, thisdist);

//				J[counti] = cal_J(J_M, tempzx, tempzy);
//				gradx[counti] = cal_interpgrad(1, tstarx, tstary);
//				grady[counti] = cal_interpgrad(2, tstarx, tstary);
				cal_interpgrad(tstarx, tstary, DIST_DX, DIST_DY, &gradx[counti], &grady[counti]);
				gradnorm = sqrt(gradx[counti] * gradx[counti] + grady[counti]*grady[counti]);
				if (fabs(gradnorm - 1.0) > 0.1)
				{
					printf("At (%lf, %lf), grad (%lf, %lf), norm %lf.\n", \
					tstarx, tstary, gradx[counti], grady[counti], gradnorm);
					exit(0);
				}

				J[counti] = 1. + thisdist * CURV[i][j];
				thiscurv = cal_curv(tstarx, tstary, DIST);
				regfactor[counti] = thiscurv/(4.*PI);
				result[counti] = HSQUARED * delta[counti] * J[counti];

				counti++;	// in the end should total up to size	//
			}
		}

		ini = 1;	//	DENOTING DATA INITIALIZED	//
	}


	starxi = zstarx[pi];
	staryi = zstary[pi];
	starxj = zstarx[pj];
	staryj = zstary[pj];

///////////////////////////////////////////////		faster kernel		//////////////////////////////

	xminusy1 = starxi-starxj;
	xminusy2 = staryi-staryj;

	thisnorm = sqrt( xminusy1*xminusy1 + xminusy2*xminusy2 );

	if (thisnorm > tau)
	{
		parameter= WAVE * thisnorm;
		commonterm = 0.25 * WAVE / thisnorm;
		BesselY1value = BesselY1(parameter);
//		partialx = commonterm * BesselY1value * xminusy1;
//		partialy =  commonterm * BesselY1value * xminusy2;
		BesselJ1value = BesselJ1(parameter);
//		partialximg = -1. * commonterm * BesselJ1value * xminusy1;
//		partialyimg = -1. * commonterm * BesselJ1value * xminusy2;

		if (DIRICHLET == 1)
			xmyinner_pd = xminusy1 * gradx[pj] + xminusy2 * grady[pj];
		else
			xmyinner_pd = xminusy1 * gradx[pi] + xminusy2 * grady[pi];

//////////////////////////////////////////////////////////////////////////////////////////////

//	if (DIRICHLET == 1)
//	{
//		partialx = cal_partial(1, starxi, staryi, starxj, staryj);// x coordinate real	//
//		partialy = cal_partial(2, starxi, staryi, starxj, staryj);// y coordinate real	//
//		partialximg = cal_partial(3, starxi, staryi, starxj, staryj);// x coordinate imaginary //
//		partialyimg = cal_partial(4, starxi, staryi, starxj, staryj);// y coordinate imaginary //
//	}
//	else
//	{
//		partialx = cal_partial(1, starxj, staryj, starxi, staryi);// x coordinate real	//
//		partialy = cal_partial(2, starxj, staryj, starxi, staryi);// y coordinate real	//
//		partialximg = cal_partial(3, starxj, staryj, starxi, staryi);// x coordinate imaginary //
//		partialyimg = cal_partial(4, starxj, staryj, starxi, staryi);// y coordinate imaginary //

//		partialx = cal_partial(1, starxi, staryi, starxj, staryj);// x coordinate real	//
//		partialy = cal_partial(2, starxi, staryi, starxj, staryj);// y coordinate real	//
//		partialximg = cal_partial(3, starxi, staryi, starxj, staryj);// x coordinate imaginary //
//		partialyimg = cal_partial(4, starxi, staryi, starxj, staryj);// y coordinate imaginary //
//	}


	//	for interior problem	//
		if (INTERIOR_MODE==1)
		{
			partial = (-1.0) * commonterm * BesselY1value * xmyinner_pd;	//	outer normal
			partialimg = commonterm * BesselJ1value * xmyinner_pd;
		}
		//	for exterior problem	//
		else
		{
			partial = -1. * commonterm * BesselY1value * xmyinner_pd;	//	outer normal
			partialimg = commonterm * BesselJ1value * xmyinner_pd;
		}
	}
//
//
//	if ( ( (starxi-starxj)*(starxi-starxj) + (staryi-staryj)*(staryi-staryj) )< threshold )
//	{
//		partial_reg = regfactor;
//		partialimg_reg = regfactorimg;
//	}
	else
	{
		partial = regfactor[i];
		partialimg = 0.0;
//		partial_reg = partial;
//		partialimg_reg = partialimg;
	}

	//	calculate matrix elements	//

	tempC[0] = result[pj] * partial ;		//	( h^3*deltaj*Jj ) * partial	//
	tempC[1] = result[pj] * partialimg;
//	if (REG_MODE==1)
//	{
//		tempC[0] = result[pj] * partial_reg;
//		tempC[1] = result[pj] * partialimg_reg;
//	}

	//	C = B + 1/2 I	//
	if (pi == pj)
	{
		if (DIRICHLET == 1)
			tempC[0] += 0.5;	//(delta*cal_J(J_M, wx, wy));
		else
			tempC[0] -= 0.5;
	}

//	return 0.;
	if ((tempC[0]!=tempC[0])||(tempC[0]+1==tempC[0])||(tempC[1]!=tempC[1])||(tempC[1]+1==tempC[1]))
	{
		printf("i = %d, j = %d, partial = %lf, C = %.10lf, Cimg = %.10lf.\n", \
			pi, pj, partial, tempC[0], tempC[1]);
		exit(0);
	}
	return 0.;
}

