
#include "BICGSTAB.h"


int CompBiCG(int size, double *density, double *densityimg, double **Kernel, double **Kernelimg, double *b, double *bimg)
{
	FILE *fp3Dx, *fprealb, *fp3Dximg, *fpimgb;
	FILE *fpverifyb, *fpverifybimg;
	FILE *fp3DC;		//	DEBUG	//	
	FILE *fp3DCimg;		//	DEBUG	//

	int i, j, k, counter;
	int compsize = 2 * size;
	double nonhomogeneous = 0.0;
	double nonhomogeneousimg = 0.0;
	double zx, zy;
//	double epsilon = FACTOR * h;
	int count = 0;;

	double tempC[2];
//	double r[size], rhat[size], x[size], v[size], p[size], s[size], t[size];
//	double rimg[size], rhatimg[size], ximg[size], vimg[size], pimg[size], simg[size], timg[size];
	double *x, *ximg, *r, *rimg, *rtilde, *rtildeimg;
	double *p, *pimg, *ptilde, *ptildeimg, *pold, *poldimg, *ptildeold, *ptildeoldimg;
	double *Akp, *Akpimg, *Akptilde, *Akptildeimg;

//	double b[3] = TESTB;
//	double bimg[3] = TESTBIMG;
//	double b[size], bimg[size];

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

	if (INTERACTIVE_ACCURACY == 1)
	{
		printf("Enter the accuracy in 10^x ( 0 > x > -10)\n");
		scanf("%lf", &accpower);
		if ( (accpower > -10.0) && (accpower < 0.0) )
			accuracy = pow(10, accpower);
		printf("Using accuracy = %.8E.\n", accuracy);
	}


	x = (double *) malloc(size * sizeof(double));
	p = (double *) malloc(size * sizeof(double));
	pold = (double *) malloc(size * sizeof(double));
	r = (double *) malloc(size * sizeof(double));
	ptilde = (double *) malloc(size * sizeof(double));
	ptildeold = (double *) malloc(size * sizeof(double));
	rtilde = (double *) malloc(size * sizeof(double));
	Akp = (double *) malloc(size * sizeof(double));
	Akptilde = (double *) malloc(size * sizeof(double));
	error = (double *) malloc(size * sizeof(double));
	ximg = (double *) malloc(size * sizeof(double));
	pimg = (double *) malloc(size * sizeof(double));
	poldimg = (double *) malloc(size * sizeof(double));
	rimg = (double *) malloc(size * sizeof(double));
	ptildeimg = (double *) malloc(size * sizeof(double));
	ptildeoldimg = (double *) malloc(size * sizeof(double));
	rtildeimg = (double *) malloc(size * sizeof(double));
	Akpimg = (double *) malloc(size * sizeof(double));
	Akptildeimg = (double *) malloc(size * sizeof(double));
	errorimg = (double *) malloc(size * sizeof(double));

	//	INITIALIZING THE THINGS NEEDED FOR BICONJUGATE GRADIENT STABILIZED	//
////////////////////////////////////////////////////////////////////////////
//	TESTB;	
//	bimg[3] = TESTBIMG;

//	fprintf(fpverifyb, "%d\n", size);
//	fprintf(fpverifybimg, "%d\n", size);
//	for (int i = 0; i < size ; i++)
//	{
//		fprintf(fpverifyb, "%lf ", b[i]);
//		fprintf(fpverifybimg, "%lf ", bimg[i]);
//	}
////////////////////////////////////////////////////////////////////////////


//	normb = sqrt(inner(size, b, b) + inner(size, bimg, bimg));
	normb = sqrt(inner(size, b, b) + inner(size, bimg, bimg));
	printf("b norm is %lf\t", normb);

	for (counter = 0; counter < size; counter++)
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
//	innerrold = inner(size, r, rtilde) + inner(size, rimg, rtildeimg);
//	innerroldimg = inner(size, rimg, rtilde) - inner(size, r, rtildeimg);
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
			density[i] = 0.;
			densityimg[i] = 0.;
		}
		return 0;
	}


	for ( i = 1; i <= 200; i++)
	{
		loop_initime = time(NULL);

		//////	FIRST MATRIX MULTIPLICATION	////////////
		for (j = 0; j < size; j++)
		{
			sum = 0.0;
			sumimg = 0.0;
			for ( k = 0; k < size; k++)
			{
				sum += (Kernel[j][k] * p[k] - Kernelimg[j][k] * pimg[k]);		//	A * p	//
				sumimg += (Kernel[j][k] * pimg[k] + Kernelimg[j][k] * p[k]);	//	A * p	//
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


		//////	SECOND MATRIX MULTIPLICATION	/////////////
		for ( j = 0; j < size; j++)
		{
			sum = 0.0;
			sumimg = 0.0;
			for ( k = 0; k < size; k++)
			{
				sum += (Kernel[k][j]*ptilde[k] + Kernelimg[k][j]*ptildeimg[k]);
				sumimg += (Kernel[k][j]*ptildeimg[k] - Kernelimg[k][j]*ptilde[k]);
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
//		comptempa = omp_inner(size, r,rtilde, 16) + omp_inner(size, rimg, rtildeimg, 16);
//		comptempb = omp_inner(size, rimg,rtilde, 16) - omp_inner(size, r, rtildeimg, 16);
//		comptempc = omp_inner(size, Akp, ptilde, 16) + omp_inner(size, Akpimg, ptildeimg, 16);
//		comptempd = omp_inner(size, Akpimg, ptilde, 16) - omp_inner(size, Akp, ptildeimg, 16);
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
//			if ( (p[j] != p[j])|| (p[j] + 1.0 ==p[j]) )
//				printf("j = %d, r[j] = %lf, v[j] = %lf, beta = %lf.\n", j, r[j], v[j], beta);
		}


		//	BEGIN beta = (r^,r)_k+1/(r^, r)k		//
		comptempc = innerrold;
		comptempd = innerroldimg;	
		comptempa = inner(size, r, rtilde) + inner(size, rimg, rtildeimg);
		comptempb = inner(size, rimg, rtilde) - inner(size, r, rtildeimg);
//		comptempa = omp_inner(size, r, rtilde, 16) + omp_inner(size, rimg, rtildeimg, 16);
//		comptempb = omp_inner(size, rimg, rtilde, 16) - omp_inner(size, r, rtildeimg, 16);
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

		for (j = 0; j < size; j++)
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
//		errornorm = sqrt(omp_inner(size, r, r, 16) + omp_inner(size, rimg, rimg, 16))/normb;
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
		density[i] = x[i];
		densityimg[i] = ximg[i];
	}
	
	free(x);
	free(p);
	free(r);
	free(pold);
	free(ptilde);
	free(rtilde);
	free(ptildeold);
	free(Akp);
	free(Akptilde);
	free(error);
	free(ximg);
	free(pimg);
	free(rimg);
	free(poldimg);
	free(ptildeimg);
	free(rtildeimg);
	free(ptildeoldimg);
	free(Akpimg);
	free(Akptildeimg);
	free(errorimg);
	return 0;
}




int omp_CompBiCG(int size, double *density, double *densityimg, double **Kernel, double **Kernelimg, \
	double *b, double *bimg, int chunk)
{
	FILE *fp3Dx, *fprealb, *fp3Dximg, *fpimgb;
	FILE *fpverifyb, *fpverifybimg;
	FILE *fp3DC;		//	DEBUG	//	
	FILE *fp3DCimg;		//	DEBUG	//

	int i, j, k, counter;
	int compsize = 2 * size;
	double nonhomogeneous = 0.0;
	double nonhomogeneousimg = 0.0;
	double zx, zy;
//	double epsilon = FACTOR * h;
	int count = 0;;

	double tempC[2];
//	double r[size], rhat[size], x[size], v[size], p[size], s[size], t[size];
//	double rimg[size], rhatimg[size], ximg[size], vimg[size], pimg[size], simg[size], timg[size];

	double *x, *ximg, *r, *rimg, *rtilde, *rtildeimg;
	double *p, *pimg, *ptilde, *ptildeimg, *pold, *poldimg, *ptildeold, *ptildeoldimg;
	double *Akp, *Akpimg, *Akptilde, *Akptildeimg;

//	double b[3] = TESTB;
//	double bimg[3] = TESTBIMG;
//	double b[size], bimg[size];

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

	if (INTERACTIVE_ACCURACY == 1)
	{
		printf("Enter the accuracy in 10^x ( 0 > x > -10)\n");
		scanf("%lf", &accpower);
		if ( (accpower > -10.0) && (accpower < 0.0) )
			accuracy = pow(10, accpower);
		printf("Using accuracy = %.8E.\n", accuracy);
	}

	x = (double *) malloc(size * sizeof(double));
	p = (double *) malloc(size * sizeof(double));
	pold = (double *) malloc(size * sizeof(double));
	r = (double *) malloc(size * sizeof(double));
	ptilde = (double *) malloc(size * sizeof(double));
	ptildeold = (double *) malloc(size * sizeof(double));
	rtilde = (double *) malloc(size * sizeof(double));
	Akp = (double *) malloc(size * sizeof(double));
	Akptilde = (double *) malloc(size * sizeof(double));
	error = (double *) malloc(size * sizeof(double));
	ximg = (double *) malloc(size * sizeof(double));
	pimg = (double *) malloc(size * sizeof(double));
	poldimg = (double *) malloc(size * sizeof(double));
	rimg = (double *) malloc(size * sizeof(double));
	ptildeimg = (double *) malloc(size * sizeof(double));
	ptildeoldimg = (double *) malloc(size * sizeof(double));
	rtildeimg = (double *) malloc(size * sizeof(double));
	Akpimg = (double *) malloc(size * sizeof(double));
	Akptildeimg = (double *) malloc(size * sizeof(double));
	errorimg = (double *) malloc(size * sizeof(double));


	//	INITIALIZING THE THINGS NEEDED FOR BICONJUGATE GRADIENT STABILIZED	//
////////////////////////////////////////////////////////////////////////////
//	TESTB;	
//	bimg[3] = TESTBIMG;

//	fprintf(fpverifyb, "%d\n", size);
//	fprintf(fpverifybimg, "%d\n", size);
//	for (int i = 0; i < size ; i++)
//	{
//		fprintf(fpverifyb, "%lf ", b[i]);
//		fprintf(fpverifybimg, "%lf ", bimg[i]);
//	}
////////////////////////////////////////////////////////////////////////////

	normb = sqrt(omp_inner(size, b, b, chunk) + omp_inner(size, bimg, bimg, chunk));
	printf("b norm is %lf\t", normb);

	#pragma omp parallel for private(counter) schedule(static,chunk)
	for (counter = 0; counter < size; counter++)
	{
		x[counter] = 0.0;
		ximg[counter] = 0.0;
		r[counter] = b[counter];
		rimg[counter] = bimg[counter];
		p[counter] = b[counter];
		pimg[counter] = bimg[counter];

		rtilde[counter] = b[counter];
		rtildeimg[counter] = (-1.0)*bimg[counter];
		ptilde[counter] = b[counter];
		ptildeimg[counter] = (-1.0)*bimg[counter];

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


	#pragma omp parallel for private(i) schedule(static, chunk) reduction(+:BC_sum)
	for ( i = 0; i < size; i++)
	{
		BC_sum += fabs(r[i]) + fabs(rimg[i]);
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
			density[i] = 0.;
			densityimg[i] = 0.;
		}
		return 0;
	}


	for ( i = 1; i <= 200; i++)
	{
		loop_initime = time(NULL);

		//////	FIRST MATRIX MULTIPLICATION	////////////
		#pragma omp parallel for private(j)
		for (j = 0; j < size; j++)
		{
			Akp[j] = 0.0;
			Akptilde[j] = 0.0;
			Akpimg[j] = 0.0;
			Akptildeimg[j] = 0.0;
		}

		#pragma omp parallel for private(j,k) schedule(static,chunk)
		for ( j = 0; j < size; j++)
		{
			for ( k = 0; k < size; k++)
			{
				if ( (Kernel[j][k] != Kernel[j][k]) || (Kernelimg[j][k] != Kernelimg[j][k]) )
				{
					printf("j, k= %d, %d, Kernel = %lf +  %lfi\n", \
					j, k, Kernel[j][k], Kernelimg[j][k]);
					exit(0);
				}
				Akp[j] += (Kernel[j][k] * p[k] - Kernelimg[j][k] * pimg[k]);
				Akpimg[j] += (Kernel[j][k] * pimg[k] + Kernelimg[j][k] * p[k]);
			}
//			if ((Akp[j] != Akp[j])||(Akp[j] +1 == Akp[j])||(Akpimg[j] != Akpimg[j])||(Akpimg[j]+1==Akpimg[j]))
//			{
//				printf("Akp wrong, j = %d, Akp[j] = %lf, Akpimg[j] = %lf.\n", j, Akp[j], Akpimg[j]);
//				exit(0);
//			}
//			printf("v[%d] = %.15lf\n", j, v[j]);
		}

//		loop_endtime = time(NULL);
//		seconds = loop_endtime - loop_initime;

//		if (i == 1)
//			printf("Made to 1st. Took %ld minutes %ld seconds.\n", seconds/60, seconds % 60);
		/////	FIRST MATRIX MULTIPLICATION ENDS	////


		//////	SECOND MATRIX MULTIPLICATION	/////////////

		#pragma omp parallel for private(j,k) schedule(static, chunk)
		for ( j = 0; j < size; j++)
		{
			for ( k = 0; k < size; k++)
			{
				Akptilde[j] += (Kernel[k][j]*ptilde[k] + Kernelimg[k][j]*ptildeimg[k]);
				Akptildeimg[j] += (Kernel[k][j]*ptildeimg[k] - Kernelimg[k][j]*ptilde[k]);
			}
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

//		#pragma omp parallel sections
//		{
//			#pragma omp section
//			{
		comptempa = omp_inner(size, r, rtilde, chunk) + omp_inner(size, rimg, rtildeimg, chunk);
		comptempb = omp_inner(size, rimg,rtilde, chunk) - omp_inner(size, r, rtildeimg, chunk);
//			}
//			#pragma omp section
//			{
		comptempc = omp_inner(size, Akp, ptilde, chunk) + omp_inner(size, Akpimg, ptildeimg, chunk);
		comptempd = omp_inner(size, Akpimg, ptilde, chunk) - omp_inner(size, Akp, ptildeimg, chunk);
//			}
//		}

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
				for (j = 0; j < size; j++)
					rtilde[j] = rtilde[j] - alpha * Akptilde[j] - alphaimg * Akptildeimg[j];
			}
			#pragma omp section
			{
				for (j = 0; j < size; j++)
					rtildeimg[j] = rtildeimg[j] + alphaimg * Akptilde[j] - alpha * Akptildeimg[j];
			}
//			if ( (p[j] != p[j])|| (p[j] + 1.0 ==p[j]) )
//				printf("j = %d, r[j] = %lf, v[j] = %lf, beta = %lf.\n", j, r[j], v[j], beta);
		}


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
		for (j = 0; j < size; j++)
		{
			pold[j] = p[j];
			poldimg[j] = pimg[j];
			ptildeold[j] = ptilde[j];
			ptildeoldimg[j] = ptildeimg[j];
		}
//
//		#pragma omp parallel sections private(j)
//		{
//			#pragma omp section
//			{
//				for ( j = 0; j < size; j++)
//					pold[j] = p[j];
//			}
//			#pragma omp section
//			{
//				for ( j = 0; j < size; j++)
//					poldimg[j] = pimg[j];
//			}
//			#pragma omp section
//			{
//				for ( j = 0; j < size; j++)
//					ptildeold[j] = ptilde[j];
//			}
//			#pragma omp section
//			{
//				for ( j = 0; j < size; j++)
//					ptildeoldimg[j] = ptildeimg[j];
//			}
//		}

		#pragma omp parallel for private(j) schedule(static, chunk)
		for (j = 0; j < size; j++)
		{
			p[j] = r[j] + beta * pold[j] - betaimg * poldimg[j];//	p = r + beta * p	//
			pimg[j] = rimg[j] + betaimg * pold[j] + beta * poldimg[j];
			ptilde[j] = rtilde[j] + beta * ptildeold[j] + betaimg * ptildeoldimg[j];
			ptildeimg[j] = rtildeimg[j] - betaimg * ptildeold[j] + beta * ptildeoldimg[j];
//			if (( p[j] != p[j]) || (p[j]==p[j]+1) )
//			{
//				printf("s wrong. j = %d, alpha = %lf, r = %lf, rimg = %lf, p = %lf\n", \
//				j, alpha, r[j], rimg[j], p[j]);
//				exit(0);
//			}
		}


//		#pragma omp parallel sections private(j)
//		{
//			#pragma omp section
//			{
//				for (j = 0; j < size; j++)
//					p[j] = r[j] + beta * pold[j] - betaimg * poldimg[j];//	p = r + beta * p	//
//			}
//			#pragma omp section
//			{
//				for (j = 0; j < size; j++)
//					pimg[j] = rimg[j] + betaimg * pold[j] + beta * poldimg[j];
//			}
//			#pragma omp section
//			{
//				for (j = 0; j < size; j++)
//					ptilde[j] = rtilde[j] + beta * ptildeold[j] + betaimg * ptildeoldimg[j];
//			}
//			#pragma omp section
//			{
//				for (j = 0; j < size; j++)
//					ptildeimg[j] = rtildeimg[j] - betaimg * ptildeold[j] + beta * ptildeoldimg[j];
//			}
//		}

//		for (j = 0; j < size; j++)
//		{
//			p[j] = r[j] + beta * pold[j] - betaimg * poldimg[j];//	p = r + beta * p	//
//			pimg[j] = rimg[j] + betaimg * pold[j] + beta * poldimg[j];
//			ptilde[j] = rtilde[j] + beta * ptildeold[j] + betaimg * ptildeoldimg[j];
//			ptildeimg[j] = rtildeimg[j] - betaimg * ptildeold[j] + beta * ptildeoldimg[j];
//		}



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
	
	#pragma omp parallel for private(i)
	for ( i = 0; i < size; i++)
	{
		density[i] = x[i];
		densityimg[i] = ximg[i];
	}

	free(x);
	free(p);
	free(r);
	free(pold);
	free(ptilde);
	free(rtilde);
	free(ptildeold);
	free(Akp);
	free(Akptilde);
	free(error);
	free(ximg);
	free(pimg);
	free(rimg);
	free(poldimg);
	free(ptildeimg);
	free(rtildeimg);
	free(ptildeoldimg);
	free(Akpimg);
	free(Akptildeimg);
	free(errorimg);	
	return 0;
}
