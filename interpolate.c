
#include "interpolate.h"


double cal_interpgrad(double xj, double yj, double **gradx, double **grady, double *interpolatex, double *interpolatey)
{
	static double x, y, lowx, lowy;
	static double temphhx, temphlx, templhx, templlx, temphhy, temphly, templhy, templly, temp, temp2;
	static double errll, errlh, errhl, errhh;
	static int i, j;

	static int flag = 0;
	

	i = (int)((yj-INT_L)/H);
	j = (int)((xj-INT_L)/H);

	lowx = xj - INT_L - (double)(j) * H;
	lowy = yj - INT_L - (double)(i) * H;
	
//	if (flag < 0)
//	{
//		printf("xj = %.10lf, yj = %.10lf, x = %lf, y = %lf.\n", xj, yj, x, y);
//		printf("templl = %lf, templh = %lf, temphl = %lf, temphh = %lf.\n", templl, templh, temphl, temphh);
//		printf("gradient of d = (%lf, %lf).\n", interpolatex, interpolatey);
//		flag++;
//	}

	templlx = gradx[i][j];
	templhx = gradx[i+1][j];
	temphlx = gradx[i][j+1];
	temphhx = gradx[i+1][j+1];


	templly = grady[i][j];
	templhy = grady[i+1][j];
	temphly = grady[i][j+1];
	temphhy = grady[i+1][j+1];


	if ( ( (templlx > 0) && (templhx > 0) && (temphlx > 0) && (temphhx > 0)) || \
		((templlx<0) && (temphlx < 0) && (templhx < 0) && (temphhx < 0)) )
	{
		*interpolatex = ( templlx * (H-lowx)*(H-lowy) + templhx * (H-lowx) * lowy + \
					temphlx * lowx * (H-lowy) + temphhx * lowx * lowy)/(H*H);
		*interpolatey = ( templly * (H-lowx)*(H-lowy) + templhy * (H-lowx) * lowy + \
					temphly * lowx * (H-lowy) + temphhy * lowx * lowy)/(H*H);
		return 0;
	}
/*	else
	{
		errll = fabs(templlx*templlx + templly*templly - 1.0);
		errlh = fabs(templhx*templhx + templhy*templhy - 1.0);
		errhl = fabs(temphlx*temphlx + temphly*temphly - 1.0);
		errhh = fabs(temphhx*temphhx + temphhy*temphhy - 1.0);
		
		if (errll < errlh)
		{
			if (errhl < errhh)
			{
				if (errll < errhl)
				{
					*interpolatex = templlx;
					*interpolatey = templly;
					return 1;
				}
				else
				{
					*interpolatex = temphlx;
					*interpolatey = temphly;
					return 2;
				}
			}
			else
			{
				if (errll < errhh)
				{
					*interpolatex = templlx;
					*interpolatey = templly;
					return 1;
				}
				else
				{
					*interpolatex = temphhx;
					*interpolatey = temphhy;
					return 3;
				}
			}
		}
		else
		{
			if (errhl < errhh)
			{
				if (errlh < errhl)
				{
					*interpolatex = templhx;
					*interpolatey = templhy;
					return 4;
				}
				else
				{
					*interpolatex = temphlx;
					*interpolatey = temphly;
					return 2;
				}
			}
			else
			{
				if (errlh < errhh)
				{
					*interpolatex = templhx;
					*interpolatey = templhy;
					return 4;
				}
				else
				{
					*interpolatex = temphhx;
					*interpolatey = temphhy;
					return 3;
				}
			}
		}
		
	}
*/
	else if (lowx < H-lowx)
	{
		if (lowy < H-lowy)
		{
			*interpolatex = templlx;
			*interpolatey = templly;
			return 1;
		}
		else
		{
			*interpolatex = templhx;
			*interpolatey = templhy;
			return 4;
		}
		
	}
	else
	{
		if (lowy < H-lowy)
		{
			*interpolatex = temphlx;
			*interpolatey = temphly;
			return 2;
		}
		else
		{
			*interpolatex = temphhx;
			*interpolatey = temphhy;
			return 3;
		}
	}			

	printf("Weird it would get here. In cal_interpgrad() in interpolate.c\n");
	exit(0);

//////////////////////////////////////////////////////////////////////////////////

	templly = grady[i][j];
	templhy = grady[i+1][j];
	temphly = grady[i][j+1];
	temphhy = grady[i+1][j+1];

	if ( ( (templly > 0) && (templhy > 0) && (temphly > 0) && (temphhy > 0)) || \
		((templly<0) && (temphly < 0) && (templhy < 0) && (temphhy < 0)) )
		*interpolatey = ( templly * (H-lowx)*(H-lowy) + templhy * (H-lowx) * lowy + \
					temphly * lowx * (H-lowy) + temphhy * lowx * lowy)/(H*H);
	else if (lowx < H-lowx)
	{
		if (lowy < H-lowy)
			*interpolatey = templly;
		else
			*interpolatey = templhy;
	}
	else
	{
		if (lowy < H-lowy)
			*interpolatey = temphly;
		else
			*interpolatey = temphhy;
	}

//	interpolatex = (ROOT[i][j].graddx * (H-lowx)*(H-lowy) + ROOT[i+1][j].graddx * (H-lowx) * lowy + \
//			ROOT[i][j+1].graddx * lowx * (H-lowy) + ROOT[i+1][j+1].graddx * lowx * lowy)/(H*H);

	if ((*interpolatex+1==*interpolatex)||(*interpolatex!=*interpolatex) || (fabs(*interpolatex) > 1.7) || \
		(*interpolatey+1==*interpolatey)||(*interpolatey!=*interpolatey) || (fabs(*interpolatey) > 1.7)  )
	{
		printf("it's gradient (%lf, %lf), xj = %lf, yj = %lf.\n", *interpolatex, *interpolatey, xj, yj);
		printf("i = %d, j = %d, lowx = %lf, lowy = %lf.\n", i, j, lowx, lowy);
//		printf("XDERll = %lf, XDERlh = %lf, XDERhl = %lf, XDERhh = %lf\n", XDER_C(i, j, 1), \
//		XDER_C(i+1, j, 1), XDER_C(i, j+1, 1), XDER_C(i+1, j+1, 1));
		printf("XDERll = %lf, XDERlh = %lf, XDERhl = %lf, XDERhh = %lf\n", gradx[i][j], \
		gradx[i+1][j], gradx[i][j+1], gradx[i+1][j+1]);
		printf("YDERll = %lf, YDERlh = %lf, YDERhl = %lf, YDERhh = %lf\n", grady[i][j], \
		grady[i+1][j], grady[i][j+1], grady[i+1][j+1]);
		exit(0);
	}

	return *interpolatex;
}


double cal_curv(double xj, double yj, double **dist)
{
	static double x, y, lowx, lowy;
	static double temphhx, temphlx, templhx, templlx, temphhy, temphly, templhy, templly, temp, temp2;
	static double interpolatex, interpolatey;
	static int i, j;

	static int flag = 0;
	

	i = (int)((yj-INT_L)/H);
	j = (int)((xj-INT_L)/H);

	lowx = xj - INT_L - (double)(j) * H;
	lowy = yj - INT_L - (double)(i) * H;
	
	templlx = XDER2(dist, i, j, FD_MODE);
	templhx = XDER2(dist, i+1, j, FD_MODE);
	temphlx = XDER2(dist, i, j+1, FD_MODE);
	temphhx = XDER2(dist, i+1, j+1, FD_MODE);

	interpolatex = ( templlx * (H-lowx)*(H-lowy) + templhx * (H-lowx) * lowy + \
			temphlx * lowx * (H-lowy) + temphhx * lowx * lowy);

	templly = YDER2(dist, i, j, FD_MODE);
	templhy = YDER2(dist, i+1, j, FD_MODE);
	temphly = YDER2(dist, i, j+1, FD_MODE);
	temphhy = YDER2(dist, i+1, j+1, FD_MODE);

	interpolatey = ( templly * (H-lowx)*(H-lowy) + templhy * (H-lowx) * lowy + \
			temphly * lowx * (H-lowy) + temphhy * lowx * lowy);

	return -1.0 * (interpolatex + interpolatey)/(H*H);

}


double cal_dist(double x, double y, double **dist, double **gradx, double **grady)
{
	static int indx, indy;
	static double lowx, lowy, smallestxx, smallestyy, tempd2, tempd1, tempd11, tempd12, tempd21, tempd22;

	indx = (int)((x - INT_L)/H);
	indy = (int)((y - INT_L)/H);

//	printf("indx = %d, indy = %d\n", indx, indy);

//	smallestxx = LARGENUM;
//	smallestyy = LARGENUM;

	lowx = x - INT_L - indx * H;
	lowy = y - INT_L - indy * H;

	tempd11 = XDER2(dist, indy, indx, FD_MODE);
	tempd21 = XDER2(dist, indy+1, indx, FD_MODE);
	tempd12 = XDER2(dist, indy, indx+1, FD_MODE);
	tempd22 = XDER2(dist, indy+1, indx+1, FD_MODE);
	if (fabs(tempd11) < fabs(tempd12))
		tempd1 = tempd11;
	else
		tempd1 = tempd12;
	if (fabs(tempd21) < fabs(tempd22) )
		tempd2 = tempd21;
	else
		tempd2 = tempd22;
	if (fabs(tempd1) < fabs(tempd2))
		smallestxx = tempd1;
	else
		smallestxx = tempd2;

	tempd11 = YDER2(dist, indy, indx, FD_MODE);
	tempd21 = YDER2(dist, indy+1, indx, FD_MODE);
	tempd12 = YDER2(dist, indy, indx+1, FD_MODE);
	tempd22 = YDER2(dist, indy+1, indx+1, FD_MODE);
	if (fabs(tempd11) < fabs(tempd12))
		tempd1 = tempd11;
	else
		tempd1 = tempd12;
	if (fabs(tempd21) < fabs(tempd22) )
		tempd2 = tempd21;
	else
		tempd2 = tempd22;
	if (fabs(tempd1) < fabs(tempd2))
		smallestyy = tempd1;
	else
		smallestyy = tempd2;

//	if (fabs(smallestxx) + fabs(smallestyy) < 5.0)
		return (dist[indy][indx] * (H-lowx) * (H-lowy) + dist[indy+1][indx] * (H-lowx) * lowy + \
			dist[indy][indx+1] * (lowx) * (H-lowy) + dist[indy+1][indx+1] * lowx * lowy)/(H*H) - \
			0.5 * ( smallestxx * (H-lowx) * lowx + smallestyy * (H-lowy) * lowy);
//	else
		return (dist[indy][indx] * (H-lowx) * (H-lowy) + dist[indy+1][indx] * (H-lowx) * lowy + \
			dist[indy][indx+1] * (lowx) * (H-lowy) + dist[indy+1][indx+1] * lowx * lowy)/(H*H);
}

