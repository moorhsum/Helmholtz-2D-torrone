
#include "fd_der.h"

double XDER_C(double **level, int i, int j, int mode)
{
	if ((j == 0)||(j == N))
		return 0.0;
	else
		return (level[i][j+1] - level[i][j-1])/(2.0*H);
}

double YDER_C(double **level, int i, int j, int mode)
{
	if ((i == 0)|| (i == N))
		return 0.0;
	else
		return (level[i+1][j] - level[i-1][j])/(2.0*H);
}



//	Computes the second derivative of the distance function in the x direction.
double XDER2(double **level, int i, int j, int mode)
{
//	static POINT *temp;
	static double DDxp, DDxn, DDxc, Dxpp, Dxp, Dxn, Dxnn;
//	temp = ROOT[i];

	if (mode == 1)
	{
		if (j < 0)
			return (level[i][abs(j)+1] - 2.0 * level[i][abs(j)] + level[i][abs(j)-1])/(H*H);
		else if (j == 0)
			return 2.0*(level[i][j+1] - level[i][j])/(H*H);
		else if (j == N)
			return 2.0*(level[i][j-1] - level[i][j])/(H*H);
		else if (j > N)
			return (level[i][2*N-j+1] - 2.0 * level[i][2*N-j] + level[i][2*N-j-1])/(H*H);
		else if (j < 3)
		{
			DDxp = (2.0*level[i][j] - 5.0*level[i][j+1] + 4.0*level[i][j+2]-level[i][j+3])/(H*H);
			DDxc = (level[i][j+1] - 2.0*level[i][j] + level[i][j-1])/(H*H);
			return (fabs(DDxp)<fabs(DDxc))?DDxp:DDxc;
		}
		else if (j > N - 3)
		{
			DDxn = (2.0*level[i][j] - 5.0*level[i][j-1] + 4.0*level[i][j-2]-level[i][j-3])/(H*H);
			DDxc = (level[i][j+1] - 2.0*level[i][j] + level[i][j-1])/(H*H);
			return (fabs(DDxn)<fabs(DDxc))?DDxn:DDxc;
		}
		else
		{
			Dxpp = level[i][j+2] - level[i][j+1];
			Dxp = level[i][j+1] - level[i][j];
			Dxn = level[i][j] - level[i][j-1];
			Dxnn = level[i][j-1] - level[i][j-2];
//			Dyp = YDER_P(i, j, 1, 0);
//			Dyn = YDER_M(i, j, 1, 0);

			DDxp = (2.0*level[i][j] - 5.0*level[i][j+1] + 4.0*level[i][j+2] - level[i][j+3])/(H*H);
			DDxc = (level[i][j+1] - 2.0*level[i][j] + level[i][j-1])/(H*H);
			DDxn = (2.0*level[i][j] - 5.0*level[i][j-1]+4.0*level[i][j-2] - level[i][j-3])/(H*H);

			if ( (Dxp * Dxnn <= 0.0) && (Dxpp * Dxn <= 0.0) )
				return (fabs(Dxp)>fabs(Dxn))?DDxp:DDxn;
			else if (Dxp * Dxnn <= 0.0)
				return DDxp;
			else if (Dxpp * Dxn <= 0.0)
				return DDxn;
			else
				return DDxc;


//			if ( (Dxp * Dxnn <= 0.0) || (Dxpp * Dxn <= 0.0) )
//			{
//				return min_abs(DDxp, DDxc, DDxn);
//				return (fabs(Dxp)<fabs(Dxn))?DDxn:DDxp;
//			}
//			else
//				return DDxc;
		}
	}
	else
	{
		if (j < 0)
			return (level[i][abs(j)+1] - 2.0 * level[i][abs(j)] + level[i][abs(j)-1])/(H*H);
		if (j == 0)
			return 2.0*(level[i][j+1] - level[i][j])/(H*H);
		else if (j == N)
			return 2.0*(level[i][j-1] - level[i][j])/(H*H);
		else if (j > N)
			return (level[i][2*N-j+1] - 2.0 * level[i][2*N-j] + level[i][2*N-j-1])/(H*H);
		else
			return (level[i][j+1] - 2.0*level[i][j] + level[i][j-1])/(H*H);
	}
}


double YDER2(double **level, int i, int j, int mode)
{
	static double DDyp, DDyc, DDyn;
	static double Dypp, Dyp, Dyn, Dynn;
	if (mode == 1)
	{
		if (i < 0)	//	Periodic condition, mirror effect	//
			return (level[abs(i)+1][j] - 2.0 * level[abs(i)][j] + level[abs(i)-1][j])/(H*H);
		else if (i == 0)
			return 2.0*( level[i+1][j] - level[i][j])/(H*H);
		else if (i == N)
			return 2.0*( level[i-1][j] - level[i][j])/(H*H);
		else if (i > N)
			return ( level[2*N-i+1][j] - 2.0 * level[2*N-i][j] + level[2*N-i-1][j])/(H*H);
		else if (i < 3)
		{
			DDyp = (2.0*level[i][j] - 5.0*level[i+1][j] + 4.0*level[i+2][j] - level[i+3][j])/(H*H);
			DDyc = (level[i+1][j] - 2.0*level[i][j] + level[i-1][j])/(H*H);
			return (fabs(DDyp)<fabs(DDyc))?DDyp:DDyc;
		}
		else if (i > N - 3)
		{
			DDyn = (2.0*level[i][j] - 5.0*level[i-1][j] + 4.0*level[i-2][j] - level[i-3][j])/(H*H);
			DDyc = (level[i+1][j] - 2.0*level[i][j] + level[i-1][j])/(H*H);
			return (fabs(DDyn)<fabs(DDyc))?DDyn:DDyc;
		}
		else
		{
			Dypp = level[i+2][j] - level[i+1][j];
			Dyp = level[i+1][j] - level[i][j];
			Dyn = level[i][j] - level[i-1][j];
			Dynn = level[i-1][j] - level[i-2][j];

			DDyp = (2.0*level[i][j] - 5.0*level[i+1][j] + 4.0*level[i+2][j] - level[i+3][j])/(H*H);
			DDyc = (level[i+1][j] - 2.0*level[i][j] + level[i-1][j])/(H*H);
			DDyn = (2.0*level[i][j] - 5.0*level[i-1][j] + 4.0*level[i-2][j] - level[i-3][j])/(H*H);

			if ( (Dyp * Dynn <= 0.0) && (Dypp * Dyn <= 0.0) )
				return (fabs(Dyp)>fabs(Dyn))?DDyp:DDyn;
			else if (Dyp * Dynn <= 0.0)
				return DDyp;
			else if (Dypp * Dyn <= 0.0)
				return DDyn;
			else
				return DDyc;



//			if ( (Dypp * Dyn <= 0.0) || (Dyp * Dynn <= 0.0) )
//			{
//				return min_abs(DDyp, DDyc, DDyn);
//				return (fabs(Dyp)>fabs(Dyn))?DDyp:DDyn;
//			}
//			else
//				return DDyc;
		}
	}
	else
	{
		if (i < 0)	//	Periodic condition	//
			return (level[abs(i)+1][j] - 2.0 * level[abs(i)][j] + level[abs(i)-1][j])/(H*H);
		else if (i == 0)
			return 2.0*(level[i+1][j] - level[i][j])/(H*H);
		else if (i == N)
			return 2.0*(level[i-1][j] - level[i][j])/(H*H);
		else if (i > N)
			return (level[2*N-i+1][j] - 2.0 * level[2*N-i][j] + level[2*N-i-1][j])/(H*H);
		else
			return (level[i+1][j] - 2.0*level[i][j] + level[i-1][j])/(H*H);
	}
}


void FD_DX(double **phi, int i, int j, double *Dxp, double *Dxn, int mode)
{
//	double **phii = phi[i];
	if (mode >= 2)
		mode = 2;
	if (mode == 1)
	{
		if (j == 0)
		{
			*Dxp = (phi[i][j+1] - phi[i][j])/H;
			*Dxn = -1.0 * (*Dxp);
		}
		else if (j == N)
		{
			*Dxn = (phi[i][j] - phi[i][j-1])/H;
			*Dxp = -1.0 * (*Dxn);
		}
		else
		{
			*Dxp = (phi[i][j+1] - phi[i][j])/H;
			*Dxn = (phi[i][j] - phi[i][j-1])/H;
		}
	}
	else if (mode == 2)
	{
		double DDxp, DDxc, DDxn;
		if (j == 0)
		{
			*Dxp = (phi[i][j+1] - phi[i][j])/H;
			*Dxn = -1.0 * (*Dxp);
		}		
		else if (j == N)
		{
			*Dxn = (phi[i][j] - phi[i][j-1])/H;
			*Dxp = -1.0 * (*Dxn);
		}
		else if ((j == 1) || (j == N-1))
		{
			*Dxp = (phi[i][j+1] - phi[i][j])/H;
			*Dxn = (phi[i][j] - phi[i][j-1])/H;
		}
		else
		{
			DDxp = XDER2(phi, i, j+1, 2);
			DDxc = XDER2(phi, i, j, 2);
			DDxn = XDER2(phi, i, j-1, 2);
			*Dxp = (phi[i][j+1] - phi[i][j])/H - 0.5 * H * minmod(DDxp, DDxc);
			*Dxn = (phi[i][j] - phi[i][j-1])/H + 0.5 * H * minmod(DDxn, DDxc);
		}
	}
//	if (mode )

}



void FD_DY(double **phi, int i, int j, double *Dyp, double *Dyn, int mode)
{
	if (mode >= 2)
		mode = 2;
	if (mode == 1)
	{
		if (i == 0)
		{
			*Dyp = (phi[i+1][j] - phi[i][j])/H;
			*Dyn = -1.0 * (*Dyp);
		}
		else if (i == N)
		{
			*Dyn = (phi[i][j] - phi[i-1][j])/H;
			*Dyp = -1.0 * (*Dyn);
		}
		else
		{
			*Dyp = (phi[i+1][j] - phi[i][j])/H;
			*Dyn = (phi[i][j] - phi[i-1][j])/H;
		}
	}
	else if (mode == 2)
	{
		double DDyp, DDyn, DDyc;
		if (i == 0)
		{
			*Dyp = (phi[i+1][j] - phi[i][j])/H;
			*Dyn = -1.0 * (*Dyp);
		}
		else if ( (i == 1)|| (i == N-1) )
		{
			*Dyp = (phi[i+1][j] - phi[i][j])/H;
			*Dyn = (phi[i][j] - phi[i-1][j])/H;
		}
		else if (i == N)
		{
			*Dyn = (phi[i][j] - phi[i-1][j])/H;
			*Dyp = -1.0 * (*Dyn);
		}
		else
		{
			DDyp = YDER2(phi, i+1, j, 2);
			DDyc = YDER2(phi, i, j, 2);
			DDyn = YDER2(phi, i-1, j, 2);
			*Dyp = (phi[i+1][j] - phi[i][j])/H - 0.5 * H * minmod(DDyp, DDyc);
			*Dyn = (phi[i][j] - phi[i-1][j])/H + 0.5 * H * minmod(DDyn, DDyc);
		}
	}
}


