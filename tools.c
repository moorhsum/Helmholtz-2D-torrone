
#include "tools.h"


double inner(int size, double vectorA[], double vectorB[])
{
	int i;
	double sum = 0.0;
	for (i = 0; i < size; i++)
	{
		sum += vectorA[i] * vectorB[i];
//		if (sum!=sum)
//		{
//			printf("i = %d, ai = %lf, bi = %lf.\n", i, vectorA[i], vectorB[i]);
//			exit(0);
//		}
	}
	return sum;
}

double omp_inner(int size, double vectorA[], double vectorB[], int chunk)
{
	int i;
	double sum = 0.0;
	#pragma omp parallel for private(i) schedule(static, chunk) reduction(+:sum)
	for (i = 0; i < size; i++)
	{
		sum += vectorA[i] * vectorB[i];
	}
	return sum;
}

double norm(int size, double vectorA[])
{	return sqrt(fabs(inner(size, vectorA, vectorA)));	}

double BesselJ0(double x)
{
	static const double polycoef[7] = J0COEFFICIENT;
	static const double fcoef[6] = F0COEFFICIENT;
	static const double thetacoef[6] = THETA0COEFFICIENT;

	int i, j;
	double sum = 0.0;

	if (x < 3.0)
	{
		if (x <= 1.0)
		{
			for (i = 0; i < 7; i++)
			{
				sum += (pow(x/3.0, 2.*i) * polycoef[i]);
			}
			return sum;
		}
		else
		{
//			double sum = 0.0;
			for (i = 0; i < 7; i++)
			{
				sum += (pow(x/3.0, 2.*i) * polycoef[i]);
			}
			return sum;
//		for (int k = 0; k < BESSELJ_APROX; k++)
//		{
//			sum += pow(-1.0, k) * pow( pow(x/2.0, k)/factorial((double)k) , 2.0);
//		}
//		return sum;
		}
	}
	else
	{
		double f = 0.0;
		double theta = 0.0;
		if (x < 5.0)
		{
			for (j = 0; j < 6; j++)
			{
				f += (pow(3.0/x, 2*j) * fcoef[j]);
				theta += (pow(3.0/x, 2*j+1) * thetacoef[j]);
			}
			theta += (x - 0.25*PI);
			return f * cos(theta)/sqrt(x);
		}
		else if (x < 12.0)
		{
			for (j = 0; j < 6; j++)
			{
				f += (pow(3.0/x, 2*j) * fcoef[j]);
				theta += (pow(3.0/x, 2*j+1) * thetacoef[j]);
			}
			theta += (x - 0.25*PI);
			return f * cos(theta)/sqrt(x);
			
		}
		else
		{
			for (j = 0; j < 6; j++)
			{
				f += (pow(3.0/x, 2*j) * fcoef[j]);
				theta += (pow(3.0/x, 2*j+1) * thetacoef[j]);
			}
			theta += (x - 0.25*PI);
			return f * cos(theta)/sqrt(x);
			
		}
	}

}

double BesselY0(double x)
{
	static const double polycoef[7] = Y0COEFFICIENT;
	static const double fcoef[6] = F0COEFFICIENT;
	static const double thetacoef[6] = THETA0COEFFICIENT;

	int i, j;
	double sum = 0.0;
//	double sum2 = EULER;
//	double sum3 = 0.0;
//	double tempJ0;
//	double result = 0.0;
//	int k = 0;
	
	if (x < 3.0)
	{
		if (x <= 0.8)
		{
			for (i = 0; i < 7; i++)
			{
				sum += (pow(x/3.0, 2*i) * polycoef[i]);
			}
			return sum + (2.0/PI) * log(x/2.0) * BesselJ0(x);
		}
		else
		{
			for (i = 0; i < 7; i++)
			{
				sum += (pow(x/3.0, 2*i) * polycoef[i]);
			}
			return sum + (2.0/PI) * log(x/2.0) * BesselJ0(x);
		}
	}
	else
	{
		double f = 0.0;
		double theta = 0.0;
		if (x < 5.0)
		{
			for (j = 0; j < 6; j++)
			{
				f += (pow(3.0/x, 2*j) * fcoef[j]);
				theta += (pow(3.0/x, 2*j+1) * thetacoef[j]);
			}
			theta += (x - 0.25*PI);
			return f * sin(theta)/sqrt(x);
		}
		else if (x < 12.0)
		{
			for (j = 0; j < 6; j++)
			{
				f += (pow(3.0/x, 2*j) * fcoef[j]);
				theta += (pow(3.0/x, 2*j+1) * thetacoef[j]);
			}
			theta += (x - 0.25*PI);
			return f * sin(theta)/sqrt(x);
		}
		else
		{
			for (j = 0; j < 6; j++)
			{
				f += (pow(3.0/x, 2*j) * fcoef[j]);
				theta += (pow(3.0/x, 2*j+1) * thetacoef[j]);
			}
			theta += (x - 0.25*PI);
			return f * sin(theta)/sqrt(x);
		}
	}

//	tempJ0 = (log( fabs(x)/2.0 ) + sum2) * BesselJ0(x);

//	for (int k = 1; k < BESSELY_APROX; k++)
//	{
//		sum3 += 1.0/k;
//		sum += pow(-1.0, k+1) * pow( pow(x/2.0, k)/factorial(k), 2.0) * sum3;
//		sum2 += 1.0/k;
//		sum += pow(-1.0, k) * pow( pow(x/2.0, k) / factorial((double) k), 2.0) * (log(fabs(x)/2.0)-sum2);
//		printf("k = %d, sum = %lf, sum2 = %lf.\n", k, sum, sum2);
//		printf("log(x/2.0) = %lf", log(fabs(x)/2.0));
//		printf("preresult = %lf, without multiply = %lf.\n", pow(pow(x/2.0,k)/factorial(k),2.0), sum + tempJ0);
//		printf("k = %d, harmonic = %lf, Current result = %lf.\n", k, sum3, (sum+tempJ0)*(2.0/PI) );
//	}
//	return sum;
//	return (sum + tempJ0)*(2.0/PI);
}

double BesselJ1(double x)
{
	static const double polycoef[7] = J1COEFFICIENT;
	static const double fcoef[6] = F1COEFFICIENT;
	static const double thetacoef[6] = THETA1COEFFICIENT;

	int i, j;
	double sum = 0.0;
	if (x <= 3.0)
	{
		for (i = 0; i < 7; i++) 
			sum += (pow(x/3.0, 2*i) * polycoef[i]);
		return sum * x;
	}
	else
	{
		double f1 = 0.0;
		double theta1 = 0.0;

		for (j = 0; j < 6; j++)
		{
			f1 += ((pow(3.0/x, 2*j) * fcoef[j]));
			theta1 += ((pow(3.0/x, 2*j + 1) * thetacoef[j]));
		}
		theta1 += (x - 0.75*PI);
		return f1 * cos(theta1)/sqrt(x);
	}
}

double BesselY1(double x)
{
	static const double polycoef[7] = Y1COEFFICIENT;
	static const double fcoef[6] = F1COEFFICIENT;
	static const double thetacoef[6] = THETA1COEFFICIENT;
	double sum = 0.0;
	int i, j;

	if (x < 3.0)
	{
		sum = (2.0/PI)*(log(x/2.0)*BesselJ1(x)-1.0/x);
		for (i = 0; i < 7; i++)
			sum += (pow(x/3.0, 2*i+1) * polycoef[i]);
		return sum;
	}
	else
	{
		double f1 = 0.0;
		double theta1 = 0.0;
		for (j = 0; j < 6; j++)
		{
			f1 += (pow(3.0/x, 2*j) * fcoef[j]);
			theta1 += (pow(3.0/x, 2*j+1) * thetacoef[j]);
		}
		theta1 += (x - 0.75*PI);
		return f1 * sin(theta1)/sqrt(x);
	}
}


double factorial(double n)
{
	if (n <= 1.0)
		return 1.0;
	else
		return n * factorial(n-1.0);
}

double max(double a, double b)
{
	if ( a >= b )
		return a;
	else
		return b;
}

double min(double a, double b)
{
	if ( a <= b )
		return a;
	else
		return b;
}

double minmod(double x, double y)
{
	if (x*y < 0)
		return 0.0;
	else
		return (fabs(x)<fabs(y))?x:y;
}

double min_abs(double a, double b, double c)
{
	double p;
	p = (fabs(a)<fabs(b))?a:b;
	return (fabs(p)<fabs(c))?p:c;
}

double Godunov(double signage, double a, double b, double c, double d)
{
	double aplus, aminus, bplus, bminus, cplus, cminus, dplus, dminus;
	aplus = max(a, 0.0) * max(a,0.0);
	aminus = min(a, 0.0) * min(a,0.0);
	bplus = max(b, 0.0) * max(b,0.0);
	bminus = min(b, 0.0) * min(b,0.0);
	cplus = max(c, 0.0) * max(c,0.0);
	cminus = min(c, 0.0) * min(c,0.0);
	dplus = max(d, 0.0) * max(d,0.0);
	dminus = min(d, 0.0) * min(d,0.0);

	if (signage > 0)
		return sqrt( max(aminus, bplus) + max(cminus, dplus));
	else if (signage < 0)
		return sqrt( max(aplus, bminus) + max(cplus, dminus));
	else
		return 0.0;
	printf("Weird Godunov\n");
	exit(0);
}
