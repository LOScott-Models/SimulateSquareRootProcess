#pragma once
#include "math.h"

void setchi2params(double dof_input, double& vb, double& vd, double& vp, int& nvexp,
	bool& is_int, bool& is_half_int)
{
//	function to set fixed parameters that are used in the repeated simulation of 
//	the chi squared distribution for cases where dof = dof -1
//	The code below is based on simulation of a gamma variate with shape parameter
//	alpha = dof/2 and beta = 1
//	Simulation method depends on the degrees of freedom parameter, alpha = dof/2

	double dof1, dof;
	int n;

	dof = dof_input - 1.0;
	vb = 0.0;
	vd = 0.0;
	vp = 0.0;
	nvexp = 0;
	is_int = 0;
	is_half_int = 0;

	if (dof <= 12.0) {
		//	Use Wallace rejection method with 2 exponential simulations
		//	Exponentials for dof = 2, 4, 6, 8, ...
		//	Gamma alpha, a = 1, 2, 3, 4, 5				
		//  Up to 6 simulations with efficiency of acceptance/rejection
		//  See Wallace, N.D. (1974), “Computer Generation of Gamma Random Variates
		//	with Non-integral Shape Parameters,” Communications of the Association
		//	for Computing Machinery (ACM) 17 (12), December 1974
		double a = dof / 2;
		dof1 = floor(a);
		n = int(dof1);
		if (fabs(dof1 - a) < 1.0e-04) is_int = 1;
		if (fabs(dof1 - a - 0.5) < 1.0e-04) is_half_int = 1;
		if (dof > 1) nvexp = n;
		dof1 = 2 * dof1;
		vp = 1.0 - 0.5 * (dof - dof1);
	}	//	end of if (dof <= 12)
	else {
		if (dof <= 25) {
			//	Cheng-Feast (1979) GKM1 ratio of uniform random numbers with acceptance/rejection
			//  See Cheng, R.C.H., and G. M. Feast (1979), “Some Simple Gamma Variate Generators,”
			//  Journal of the Royal Statistical Society, Series C (Applied Statistics), Vol. 28,
			//  No. 3 (1979):  290-295.
			//	consider GMK2 method in Cheng  Feast
			dof1 = 0.5 * dof - 1.0;
			vb = (0.5 * dof - (1.0 / (3.0 * dof))) / dof1;
			vp = 2.0 / dof1;
			vd = vp + 2.0;
		}		//		end of else for if (dof <= 20)
		else {
			//	use normal approximation for cube root(chi-squared/dof), Hilferty-Wilson method
			//	Wilson, E. B., and M. M. Hilferty. 1931. “The Distribution of Chi-Square.”
			//	Proceedings of the National Academy of Sciences. USA 17 (12): 684–688
			vd = 1.0 - 2.0 / (9.0 * dof) + 80.0 / (2187.0 * dof * dof * dof);
			vb = sqrt(2.0 / (9.0 * dof) - 104.0 / (2187.0 * dof * dof * dof));
		}	//	end of else for if dof < 25
	}		//	end of else for if dof = 12

}


unsigned int** matrix_unsigned_int(int nrow, int ncol)
// allocate an unsigned int matrix using malloc
{
	int i;
	unsigned int** m;
	// allocate pointers to rows 
	m = (unsigned int**)malloc((size_t)(nrow * sizeof(unsigned int*)));
	// allocate the necessary memory and set pointers to the rows 
	m[0] = (unsigned int*)malloc((size_t)(nrow * ncol * sizeof(unsigned int)));
	for (i = 1; i < nrow; i++) m[i] = m[i - 1] + ncol;
	// return pointer to array of pointers to rows
	return m;
}

void free_matrix_unsigned_int(unsigned int** m)
// free an unsigned int matrix allocated by matrix_unsigned_int() 
{
	free(m[0]);
	free(m);
}


//	The following funcrtion implements the step ahead algorithm for MRG32k3a
//	in Bradley, T., J. du Toit, M. Giles, R. Tong, and P. Woodhams, 
//	“Parallelisation Techniques for Random Number Generators," available
//	at https://www.nag.com/IndustryArticles/gpu_gems_article.pdf
//	This function uses divdie and conquer
void SkipAhead_MRG32k3a(int n, unsigned int An1[3][3], unsigned int An2[3][3])
{
	const unsigned int im1 = 4294967087;
	const unsigned int im2 = 4294944443;
	const unsigned int ia12 = 1403580;
	const unsigned int ia13n = 810728;
	const unsigned int ia21 = 527612;
	const unsigned int ia23n = 1370589;
	int i, j, k, ii;
	long long lp1, lp2;
	long long A1[3][3], A2[3][3], B1[3][3], B2[3][3], C1[3][3], C2[3][3];
	unsigned long long BB1[3][3], BB2[3][3], CC1[3][3], CC2[3][3];

	A1[0][0] = 0; A1[0][1] = ia12;
	A1[0][2] = 0;
	A1[0][2] -= ia13n;
	//	A1[0][2] = -ia13n;
	A1[1][0] = 1; A1[1][1] = 0; A1[1][2] = 0;
	A1[2][0] = 0; A1[2][1] = 1; A1[2][2] = 0;

	A2[0][0] = ia21; A2[0][1] = 0;
	A2[0][2] = 0;
	A2[0][2] -= ia23n;
	//	A2[0][2] = -ia23n;
	A2[1][0] = 1; A2[1][1] = 0; A2[1][2] = 0;
	A2[2][0] = 0; A2[2][1] = 1; A2[2][2] = 0;

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			B1[i][j] = A1[i][j];
			B2[i][j] = A2[i][j];
		}
	}

	for (ii = 2; ii <= 4; ii++) {
		//	pre-multiply by Ai, calculating with 64 bit signed integers
		C1[0][0] = A1[0][0] * B1[0][0] + A1[0][1] * B1[1][0] + A1[0][2] * B1[2][0];
		C1[0][1] = A1[0][0] * B1[0][1] + A1[0][1] * B1[1][1] + A1[0][2] * B1[2][1];
		C1[0][2] = A1[0][0] * B1[0][2] + A1[0][1] * B1[1][2] + A1[0][2] * B1[2][2];
		C1[1][0] = A1[1][0] * B1[0][0] + A1[1][1] * B1[1][0] + A1[1][2] * B1[2][0];
		C1[1][1] = A1[1][0] * B1[0][1] + A1[1][1] * B1[1][1] + A1[1][2] * B1[2][1];
		C1[1][2] = A1[1][0] * B1[0][2] + A1[1][1] * B1[1][2] + A1[1][2] * B1[2][2];
		C1[2][0] = A1[2][0] * B1[0][0] + A1[2][1] * B1[1][0] + A1[2][2] * B1[2][0];
		C1[2][1] = A1[2][0] * B1[0][1] + A1[2][1] * B1[1][1] + A1[2][2] * B1[2][1];
		C1[2][2] = A1[2][0] * B1[0][2] + A1[2][1] * B1[1][2] + A1[2][2] * B1[2][2];

		C2[0][0] = A2[0][0] * B2[0][0] + A2[0][1] * B2[1][0] + A2[0][2] * B2[2][0];
		C2[0][1] = A2[0][0] * B2[0][1] + A2[0][1] * B2[1][1] + A2[0][2] * B2[2][1];
		C2[0][2] = A2[0][0] * B2[0][2] + A2[0][1] * B2[1][2] + A2[0][2] * B2[2][2];
		C2[1][0] = A2[1][0] * B2[0][0] + A2[1][1] * B2[1][0] + A2[1][2] * B2[2][0];
		C2[1][1] = A2[1][0] * B2[0][1] + A2[1][1] * B2[1][1] + A2[1][2] * B2[2][1];
		C2[1][2] = A2[1][0] * B2[0][2] + A2[1][1] * B2[1][2] + A2[1][2] * B2[2][2];
		C2[2][0] = A2[2][0] * B2[0][0] + A2[2][1] * B2[1][0] + A2[2][2] * B2[2][0];
		C2[2][1] = A2[2][0] * B2[0][1] + A2[2][1] * B2[1][1] + A2[2][2] * B2[2][1];
		C2[2][2] = A2[2][0] * B2[0][2] + A2[2][1] * B2[1][2] + A2[2][2] * B2[2][2];

		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				lp1 = C1[i][j];
				lp1 = lp1 % im1;
				if (lp1 < 0) lp1 += im1;
				B1[i][j] = lp1;
				lp2 = C2[i][j];
				lp2 = lp2 % im2;
				if (lp2 < 0) lp2 += im2;
				B2[i][j] = lp2;
			}
		}

	}

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			BB1[i][j] = B1[i][j];
			BB2[i][j] = B2[i][j];
		}
	}

	ii = 8;
	while (ii <= n) {
		//	here we are squaring the matrix at each round
		CC1[0][0] = ((BB1[0][0] * BB1[0][0]) % im1 + (BB1[0][1] * BB1[1][0]) % im1 + (BB1[0][2] * BB1[2][0]) % im1) % im1;
		CC1[0][1] = ((BB1[0][0] * BB1[0][1]) % im1 + (BB1[0][1] * BB1[1][1]) % im1 + (BB1[0][2] * BB1[2][1]) % im1) % im1;
		CC1[0][2] = ((BB1[0][0] * BB1[0][2]) % im1 + (BB1[0][1] * BB1[1][2]) % im1 + (BB1[0][2] * BB1[2][2]) % im1) % im1;
		CC1[1][0] = ((BB1[1][0] * BB1[0][0]) % im1 + (BB1[1][1] * BB1[1][0]) % im1 + (BB1[1][2] * BB1[2][0]) % im1) % im1;
		CC1[1][1] = ((BB1[1][0] * BB1[0][1]) % im1 + (BB1[1][1] * BB1[1][1]) % im1 + (BB1[1][2] * BB1[2][1]) % im1) % im1;
		CC1[1][2] = ((BB1[1][0] * BB1[0][2]) % im1 + (BB1[1][1] * BB1[1][2]) % im1 + (BB1[1][2] * BB1[2][2]) % im1) % im1;
		CC1[2][0] = ((BB1[2][0] * BB1[0][0]) % im1 + (BB1[2][1] * BB1[1][0]) % im1 + (BB1[2][2] * BB1[2][0]) % im1) % im1;
		CC1[2][1] = ((BB1[2][0] * BB1[0][1]) % im1 + (BB1[2][1] * BB1[1][1]) % im1 + (BB1[2][2] * BB1[2][1]) % im1) % im1;
		CC1[2][2] = ((BB1[2][0] * BB1[0][2]) % im1 + (BB1[2][1] * BB1[1][2]) % im1 + (BB1[2][2] * BB1[2][2]) % im1) % im1;

		CC2[0][0] = ((BB2[0][0] * BB2[0][0]) % im2 + (BB2[0][1] * BB2[1][0]) % im2 + (BB2[0][2] * BB2[2][0]) % im2) % im2;
		CC2[0][1] = ((BB2[0][0] * BB2[0][1]) % im2 + (BB2[0][1] * BB2[1][1]) % im2 + (BB2[0][2] * BB2[2][1]) % im2) % im2;
		CC2[0][2] = ((BB2[0][0] * BB2[0][2]) % im2 + (BB2[0][1] * BB2[1][2]) % im2 + (BB2[0][2] * BB2[2][2]) % im2) % im2;
		CC2[1][0] = ((BB2[1][0] * BB2[0][0]) % im2 + (BB2[1][1] * BB2[1][0]) % im2 + (BB2[1][2] * BB2[2][0]) % im2) % im2;
		CC2[1][1] = ((BB2[1][0] * BB2[0][1]) % im2 + (BB2[1][1] * BB2[1][1]) % im2 + (BB2[1][2] * BB2[2][1]) % im2) % im2;
		CC2[1][2] = ((BB2[1][0] * BB2[0][2]) % im2 + (BB2[1][1] * BB2[1][2]) % im2 + (BB2[1][2] * BB2[2][2]) % im2) % im2;
		CC2[2][0] = ((BB2[2][0] * BB2[0][0]) % im2 + (BB2[2][1] * BB2[1][0]) % im2 + (BB2[2][2] * BB2[2][0]) % im2) % im2;
		CC2[2][1] = ((BB2[2][0] * BB2[0][1]) % im2 + (BB2[2][1] * BB2[1][1]) % im2 + (BB2[2][2] * BB2[2][1]) % im2) % im2;
		CC2[2][2] = ((BB2[2][0] * BB2[0][2]) % im2 + (BB2[2][1] * BB2[1][2]) % im2 + (BB2[2][2] * BB2[2][2]) % im2) % im2;

		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				BB1[i][j] = CC1[i][j];
				BB2[i][j] = CC2[i][j];
			}

		}

		ii = 2 * ii;
	}

	k = ii / 2;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			B1[i][j] = BB1[i][j];
			B2[i][j] = BB2[i][j];
		}
	}

	for (ii = (k + 1); ii <= n; ii++) {
		//	pre-multiply by Ai, calculating with 64 bit signed integers
		C1[0][0] = A1[0][0] * B1[0][0] + A1[0][1] * B1[1][0] + A1[0][2] * B1[2][0];
		C1[0][1] = A1[0][0] * B1[0][1] + A1[0][1] * B1[1][1] + A1[0][2] * B1[2][1];
		C1[0][2] = A1[0][0] * B1[0][2] + A1[0][1] * B1[1][2] + A1[0][2] * B1[2][2];
		C1[1][0] = A1[1][0] * B1[0][0] + A1[1][1] * B1[1][0] + A1[1][2] * B1[2][0];
		C1[1][1] = A1[1][0] * B1[0][1] + A1[1][1] * B1[1][1] + A1[1][2] * B1[2][1];
		C1[1][2] = A1[1][0] * B1[0][2] + A1[1][1] * B1[1][2] + A1[1][2] * B1[2][2];
		C1[2][0] = A1[2][0] * B1[0][0] + A1[2][1] * B1[1][0] + A1[2][2] * B1[2][0];
		C1[2][1] = A1[2][0] * B1[0][1] + A1[2][1] * B1[1][1] + A1[2][2] * B1[2][1];
		C1[2][2] = A1[2][0] * B1[0][2] + A1[2][1] * B1[1][2] + A1[2][2] * B1[2][2];

		C2[0][0] = A2[0][0] * B2[0][0] + A2[0][1] * B2[1][0] + A2[0][2] * B2[2][0];
		C2[0][1] = A2[0][0] * B2[0][1] + A2[0][1] * B2[1][1] + A2[0][2] * B2[2][1];
		C2[0][2] = A2[0][0] * B2[0][2] + A2[0][1] * B2[1][2] + A2[0][2] * B2[2][2];
		C2[1][0] = A2[1][0] * B2[0][0] + A2[1][1] * B2[1][0] + A2[1][2] * B2[2][0];
		C2[1][1] = A2[1][0] * B2[0][1] + A2[1][1] * B2[1][1] + A2[1][2] * B2[2][1];
		C2[1][2] = A2[1][0] * B2[0][2] + A2[1][1] * B2[1][2] + A2[1][2] * B2[2][2];
		C2[2][0] = A2[2][0] * B2[0][0] + A2[2][1] * B2[1][0] + A2[2][2] * B2[2][0];
		C2[2][1] = A2[2][0] * B2[0][1] + A2[2][1] * B2[1][1] + A2[2][2] * B2[2][1];
		C2[2][2] = A2[2][0] * B2[0][2] + A2[2][1] * B2[1][2] + A2[2][2] * B2[2][2];

		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				lp1 = C1[i][j];
				lp1 = lp1 % im1;
				if (lp1 < 0) lp1 += im1;
				B1[i][j] = lp1;
				lp2 = C2[i][j];
				lp2 = lp2 % im2;
				if (lp2 < 0) lp2 += im2;
				B2[i][j] = lp2;
			}
		}
	}

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			An1[i][j] = (unsigned int) B1[i][j];
			An2[i][j] = (unsigned int) B2[i][j];
		}
	}

	return;

}

/*
//	Following function is not used in this project
void SkipAhead2_MRG32k3a(int n, unsigned int** An1, unsigned int** An2, unsigned int** Bn1, unsigned int** Bn2)
{
	const unsigned int im1 = 4294967087;
	const unsigned int im2 = 4294944443;
	int i, j, k, ii;
	//	long long kmod, lp1, lp2;
	unsigned long long A1[3][3], A2[3][3], BB1[3][3], BB2[3][3], CC1[3][3], CC2[3][3];

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			A1[i][j] = An1[i][j];
			A2[i][j] = An2[i][j];
			BB1[i][j] = An1[i][j];
			BB2[i][j] = An2[i][j];
		}
	}

	ii = 2;
	while (ii <= n) {
		//	here we are squaring the matrix at each round
		CC1[0][0] = ((BB1[0][0] * BB1[0][0]) % im1 + (BB1[0][1] * BB1[1][0]) % im1 + (BB1[0][2] * BB1[2][0]) % im1) % im1;
		CC1[0][1] = ((BB1[0][0] * BB1[0][1]) % im1 + (BB1[0][1] * BB1[1][1]) % im1 + (BB1[0][2] * BB1[2][1]) % im1) % im1;
		CC1[0][2] = ((BB1[0][0] * BB1[0][2]) % im1 + (BB1[0][1] * BB1[1][2]) % im1 + (BB1[0][2] * BB1[2][2]) % im1) % im1;
		CC1[1][0] = ((BB1[1][0] * BB1[0][0]) % im1 + (BB1[1][1] * BB1[1][0]) % im1 + (BB1[1][2] * BB1[2][0]) % im1) % im1;
		CC1[1][1] = ((BB1[1][0] * BB1[0][1]) % im1 + (BB1[1][1] * BB1[1][1]) % im1 + (BB1[1][2] * BB1[2][1]) % im1) % im1;
		CC1[1][2] = ((BB1[1][0] * BB1[0][2]) % im1 + (BB1[1][1] * BB1[1][2]) % im1 + (BB1[1][2] * BB1[2][2]) % im1) % im1;
		CC1[2][0] = ((BB1[2][0] * BB1[0][0]) % im1 + (BB1[2][1] * BB1[1][0]) % im1 + (BB1[2][2] * BB1[2][0]) % im1) % im1;
		CC1[2][1] = ((BB1[2][0] * BB1[0][1]) % im1 + (BB1[2][1] * BB1[1][1]) % im1 + (BB1[2][2] * BB1[2][1]) % im1) % im1;
		CC1[2][2] = ((BB1[2][0] * BB1[0][2]) % im1 + (BB1[2][1] * BB1[1][2]) % im1 + (BB1[2][2] * BB1[2][2]) % im1) % im1;

		CC2[0][0] = ((BB2[0][0] * BB2[0][0]) % im2 + (BB2[0][1] * BB2[1][0]) % im2 + (BB2[0][2] * BB2[2][0]) % im2) % im2;
		CC2[0][1] = ((BB2[0][0] * BB2[0][1]) % im2 + (BB2[0][1] * BB2[1][1]) % im2 + (BB2[0][2] * BB2[2][1]) % im2) % im2;
		CC2[0][2] = ((BB2[0][0] * BB2[0][2]) % im2 + (BB2[0][1] * BB2[1][2]) % im2 + (BB2[0][2] * BB2[2][2]) % im2) % im2;
		CC2[1][0] = ((BB2[1][0] * BB2[0][0]) % im2 + (BB2[1][1] * BB2[1][0]) % im2 + (BB2[1][2] * BB2[2][0]) % im2) % im2;
		CC2[1][1] = ((BB2[1][0] * BB2[0][1]) % im2 + (BB2[1][1] * BB2[1][1]) % im2 + (BB2[1][2] * BB2[2][1]) % im2) % im2;
		CC2[1][2] = ((BB2[1][0] * BB2[0][2]) % im2 + (BB2[1][1] * BB2[1][2]) % im2 + (BB2[1][2] * BB2[2][2]) % im2) % im2;
		CC2[2][0] = ((BB2[2][0] * BB2[0][0]) % im2 + (BB2[2][1] * BB2[1][0]) % im2 + (BB2[2][2] * BB2[2][0]) % im2) % im2;
		CC2[2][1] = ((BB2[2][0] * BB2[0][1]) % im2 + (BB2[2][1] * BB2[1][1]) % im2 + (BB2[2][2] * BB2[2][1]) % im2) % im2;
		CC2[2][2] = ((BB2[2][0] * BB2[0][2]) % im2 + (BB2[2][1] * BB2[1][2]) % im2 + (BB2[2][2] * BB2[2][2]) % im2) % im2;

		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				BB1[i][j] = CC1[i][j];
				BB2[i][j] = CC2[i][j];
			}
		}

		ii = 2 * ii;
	}

	k = ii / 2;

	for (ii = (k + 1); ii <= n; ii++) {
		//	pre-multiply by Ai, calculating with 64 bit signed integers

		CC1[0][0] = ((A1[0][0] * BB1[0][0]) % im1 + (A1[0][1] * BB1[1][0]) % im1 + (A1[0][2] * BB1[2][0]) % im1) % im1;
		CC1[0][1] = ((A1[0][0] * BB1[0][1]) % im1 + (A1[0][1] * BB1[1][1]) % im1 + (A1[0][2] * BB1[2][1]) % im1) % im1;
		CC1[0][2] = ((A1[0][0] * BB1[0][2]) % im1 + (A1[0][1] * BB1[1][2]) % im1 + (A1[0][2] * BB1[2][2]) % im1) % im1;
		CC1[1][0] = ((A1[1][0] * BB1[0][0]) % im1 + (A1[1][1] * BB1[1][0]) % im1 + (A1[1][2] * BB1[2][0]) % im1) % im1;
		CC1[1][1] = ((A1[1][0] * BB1[0][1]) % im1 + (A1[1][1] * BB1[1][1]) % im1 + (A1[1][2] * BB1[2][1]) % im1) % im1;
		CC1[1][2] = ((A1[1][0] * BB1[0][2]) % im1 + (A1[1][1] * BB1[1][2]) % im1 + (A1[1][2] * BB1[2][2]) % im1) % im1;
		CC1[2][0] = ((A1[2][0] * BB1[0][0]) % im1 + (A1[2][1] * BB1[1][0]) % im1 + (A1[2][2] * BB1[2][0]) % im1) % im1;
		CC1[2][1] = ((A1[2][0] * BB1[0][1]) % im1 + (A1[2][1] * BB1[1][1]) % im1 + (A1[2][2] * BB1[2][1]) % im1) % im1;
		CC1[2][2] = ((A1[2][0] * BB1[0][2]) % im1 + (A1[2][1] * BB1[1][2]) % im1 + (A1[2][2] * BB1[2][2]) % im1) % im1;

		CC2[0][0] = ((A2[0][0] * BB2[0][0]) % im2 + (A2[0][1] * BB2[1][0]) % im2 + (A2[0][2] * BB2[2][0]) % im2) % im2;
		CC2[0][1] = ((A2[0][0] * BB2[0][1]) % im2 + (A2[0][1] * BB2[1][1]) % im2 + (A2[0][2] * BB2[2][1]) % im2) % im2;
		CC2[0][2] = ((A2[0][0] * BB2[0][2]) % im2 + (A2[0][1] * BB2[1][2]) % im2 + (A2[0][2] * BB2[2][2]) % im2) % im2;
		CC2[1][0] = ((A2[1][0] * BB2[0][0]) % im2 + (A2[1][1] * BB2[1][0]) % im2 + (A2[1][2] * BB2[2][0]) % im2) % im2;
		CC2[1][1] = ((A2[1][0] * BB2[0][1]) % im2 + (A2[1][1] * BB2[1][1]) % im2 + (A2[1][2] * BB2[2][1]) % im2) % im2;
		CC2[1][2] = ((A2[1][0] * BB2[0][2]) % im2 + (A2[1][1] * BB2[1][2]) % im2 + (A2[1][2] * BB2[2][2]) % im2) % im2;
		CC2[2][0] = ((A2[2][0] * BB2[0][0]) % im2 + (A2[2][1] * BB2[1][0]) % im2 + (A2[2][2] * BB2[2][0]) % im2) % im2;
		CC2[2][1] = ((A2[2][0] * BB2[0][1]) % im2 + (A2[2][1] * BB2[1][1]) % im2 + (A2[2][2] * BB2[2][1]) % im2) % im2;
		CC2[2][2] = ((A2[2][0] * BB2[0][2]) % im2 + (A2[2][1] * BB2[1][2]) % im2 + (A2[2][2] * BB2[2][2]) % im2) % im2;

		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				BB1[i][j] = CC1[i][j];
				BB2[i][j] = CC2[i][j];
			}
		}

	}

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			Bn1[i][j] = unsigned int(BB1[i][j]);
			Bn2[i][j] = unsigned int(BB2[i][j]);
		}
	}

	return;

}
*/
