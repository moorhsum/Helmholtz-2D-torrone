
#ifndef PDE_H_DEFINED
#define PDE_H_DEFINED


#define LAPLACE 0		//	0 for Helmholtz, 1 for Laplace		//
#define DIRICHLET 0		//	Dirichlet or Neumann condition		//
#define KITE_TEST 0		//	A nonconvex shape test, specifically for scattering problem comparison with Kress
#define INTERIOR_MODE 0		//	interior solution or exterior solution	//
#define REG_MODE 1		//	regularize or not			//
#define DIST_FILEMODE 1		//	1 to take distance function from file	//
#define DER_MODE 2		//	finite difference derivative order	//

#define POLY_TEST 0		//	test polynomial weight	//
#define SINE_TEST 1		//	test shifted sine weight	//

#define TRAD_ONESIDED 0		//	1 if KTT band only taken on one side	//
#define COMBO_ONESIDED 1	//	1 if combo approach only on one side	//
#define TAU_FACTOR 10.0 	//	tau = H/TAU_FACTOR	//

#define PRINT_KERNEL 0		//	print the kernel (respective to the sample point) or not	//
#define BOUNDARY_UTEST 0	//	1 if want to test solution accuracy on boundary	//
#define TEST_SINGLE_LAYER 1	//	Test single layer or not	//

#define PLANEWAVEX -1.0		//	the incidental plane wave direction for scattering	//
#define PLANEWAVEY 0.0		//	the incidental plane wave direction for scattering	//


#endif
