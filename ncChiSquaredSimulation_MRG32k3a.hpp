#pragma once
#include "math.h"

const double norm = 2.328306549295728e-10;
const unsigned int im1 = 4294967087;
const unsigned int im2 = 4294944443;
const unsigned int ia12 = 1403580;
const unsigned int ia13n = 810728;
const unsigned int ia21 = 527612;
const unsigned int ia23n = 1370589;

double rand_u01(unsigned int seed[]);
double sninvdev(unsigned int seed[]);
double stdNormalInverse(double p);
double rand_log_u_compl(unsigned int seed[]);
double chi2_dev(double dof, unsigned int seed[]);
double chi2dev_set(double dof, unsigned int seed[], double b, double d, double p,
	int n, bool is_int, bool is_half_int);
double rn_poisson(double u, double x, double log_Factorial[], double wgt[]);

//	Code to simulate QB method for square root process

double simulateQB(double x0, unsigned int seed[], double b, double df_term,
	double sig_dt2, double dfterm2, double c, double mu_Beta, double m2_Beta,
	double Beta_power)
{
	double x, sqrt_term, y, v, tem, tem2, quad_a, quad_b, quad_c, p_mix,
		b_mix, lambda_nc, lambda_cf;

	sqrt_term = b * x0 + df_term;
	if (df_term >= 0.0) {
		tem = sig_dt2 * sninvdev(seed);
		y = sqrt(sqrt_term);
		x = (y + tem) * (y + tem);
	}
	else {
		lambda_nc = c * b * x0;
		if (sqrt_term >= 0.0 && lambda_nc > 4.0) {
			tem = sig_dt2 * sninvdev(seed);
			y = sqrt(sqrt_term);
			x = (y + tem) * (y + tem);
		}
		else {
			//	Simulate a mixture of beta and non-central chi square(1) with moment matching
			//	mean
			y = c * (b * x0 + dfterm2);
			//	variance
			v = c * 4.0 * (b * x0 + 0.5 * dfterm2);
			lambda_cf = (lambda_nc * lambda_nc + 6.0 * lambda_nc + 3.0) /
				((lambda_nc + 1.0) * (lambda_nc + 1.0));
			quad_b = (v + y * y) + m2_Beta - 2.0 * mu_Beta * y * lambda_cf;
			quad_a = mu_Beta * mu_Beta * lambda_cf - m2_Beta;
			quad_c = y * y * (lambda_cf - 1.0) - v;
			tem = quad_b * quad_b - 4.0 * quad_a * quad_c;
			//	Check to see if this term is ever negative
			if (tem < 0.0) printf("  Warning square root term is negative %12.6e \n", tem);
			p_mix = (-quad_b + sqrt(tem)) / (2.0 * quad_a);
			if (p_mix >= 0.0 && p_mix <= 1.0) b_mix = (y - p_mix * mu_Beta) / ((1.0 - p_mix) * (lambda_nc + 1.0));
			else {
				p_mix = (-quad_b - sqrt(tem)) / (2.0 * quad_a);
				b_mix = (y - p_mix * mu_Beta) / ((1.0 - p_mix) * (lambda_nc + 1.0));
			}
			//	Check to see if a bad value for the probability is generated
			if (p_mix < 0.0 || p_mix > 1.0) printf(" warning:  p_mix %f  b_mix %f \n", p_mix, b_mix);
			//	Simulate U(0,1) for mixture
			tem = rand_u01(seed);
			if (tem <= p_mix) {
				tem2 = tem / p_mix;
				x = pow(tem2, Beta_power);
			}
			else {
				tem2 = sqrt(lambda_nc) + stdNormalInverse((1.0 - tem) / (1.0 - p_mix));
				x = b_mix * tem2 * tem2;
			}
			x = x / c;
		}	//	end of else for if (sqrt_term >= 0 && lambda_nc > 2.0)

	}	//	end of else for if (dfterm >= 0)

	return x;
}

double simulateQE(double x0, unsigned int seed[], double b, double QE_sig_dt2, 
	double dfterm2)
{	
//	Simulate a non-central chi-squared for x, using Andersen QE method
//	See Andersen, L. (2008), “Simple and efficient simulation of the Heston
//	stochastic volatility model,” Journal of Computational Finance 11 (3): 1–42
//	and Andersen, L., P. Jackel, C. Kahl, “Simulation of Square?Root Processes,”
//	chapter in Encyclopedia of Quantitative Finance, May 2010

	double x, tem, QE_m, QE_s2, QE_psi, QE_a, QE_b2, QE_b, QE_p, QE_beta;

	QE_m = b * x0 + dfterm2;
	QE_s2 = QE_sig_dt2 * (b * x0 + 0.5 * dfterm2);
	QE_psi = QE_s2 / (QE_m * QE_m);
	//	if (QE_psi <= psi_c), where psi_c = [1.0, 2.0], Andersen recommends 1.5
	if (QE_psi <= 1.5) {
		tem = 2.0 / QE_psi;
		QE_b2 = tem - 1.0 + sqrt(tem) * sqrt(tem - 1.0);
		QE_b = sqrt(QE_b2);
		QE_a = QE_m / (1.0 + QE_b2);
		tem = QE_b + sninvdev(seed);
		x = QE_a * tem * tem;
	}
	else {
		tem = 1.0 + QE_psi;
		QE_p = (QE_psi - 1.0) / tem;
		QE_beta = 2.0 / (QE_m *tem);
		tem = rand_u01(seed);
		if (tem <= QE_p) x = 0.0;
		else x = -log((1.0 - tem) / (1.0 - QE_p)) / QE_beta;
	}
	return x;
}

double nc_chi2(double x0, double dof, double c, double vrho, unsigned int seed[], 
	double b, double d, double p, int n, bool is_int, bool is_half_int, 
	double log_Factorial[], double wksp[])
{
//	This function uses exact simulation to simulate a non-central chi-squared
//	and returns  z / c for the square root process
//	See Johnson, N., S. Kotz, and N. Balakrishnan (1995), Continuous Univariate
//	Distributions, Vol. 2, 2nd Edition, Wiley Interscience
//	and Glasserman, P. (2003), Monte Carlo Methods in Financial Engineering,
//	Springer Verlag, New York

	double z, parmnc, temsq, tem, ans;
	//		Simulate a non-central chi-squared for x
	if (dof >= 0.9999) {
		parmnc = x0 * c * vrho;
		tem = sninvdev(seed);
		temsq = sqrt(parmnc);
		z = (tem + temsq) * (tem + temsq);
		if (dof > 1.0001) z = z + chi2dev_set(dof - 1.0, seed, b, d, p, n, is_int, is_half_int);
	}
	else {
		parmnc = 0.5 * x0 * c * vrho;
		if (parmnc < 200.0) {
			tem = rand_u01(seed);
			z = rn_poisson(tem, parmnc, log_Factorial, wksp);
		}
		else {
			//	for parmnc >= 200.0, use normal approximation
			tem = sninvdev(seed) * sqrt(parmnc) + parmnc;
			if (tem < 0.0) tem = 0.0;
			z = floor(tem + 0.5);
		}
		tem = dof + 2.0 * z;
		if (tem > 0.0) z = chi2_dev(tem, seed);
		else z = 0.0;
	}
	ans = z / c;
	return ans;
}

void setchi2params(double dof_input, double& vb, double& vd, double& vp, int& nvexp, 
	bool &is_int, bool &is_half_int)
{
	//	function to set chi-sqared parameters that will be reused for cases where dof = dof -1
	double dof1, dof;
	int n;

	dof = dof_input - 1.0;
	vb = 0.0;
	vd = 0.0;
//	vdof1 = 0.0;
	vp = 0.0;
	nvexp = 0; 
	is_int = 0;
	is_half_int = 0;

	if (dof <= 12.0) {
		//	Use Wallace rejection method with 2 exponential simulations
		//	Exponentials for dof = 2, 4, 6, 8, ...
		//	Gamma alpha, a = 1, 2, 3, 4, 5
		//	Compute switch probability
		double a = dof / 2;
		dof1 = floor(a);
		n = int(dof1);
		if (fabs(dof1-a) < 1.0e-04) is_int = 1;
		if(fabs(dof1-a-0.5) < 1.0e-04) is_half_int =1;
	//	vdof1 = dof;
		if (dof > 1) nvexp = n;
		dof1 = 2 * dof1;
		vp = 1.0 - 0.5 * (dof - dof1);
	}	//	end of if (dof <= 12)
	else {
		if (dof <= 25) {
			dof1 = 0.5 * dof - 1.0;
			vb = (0.5 * dof - (1.0 / (3.0 * dof))) / dof1;
			vp = 2.0 / dof1;
			vd = vp + 2.0;
		//	vdof1 = dof;
		}		//		end of else for if (dof <= 20)
		else {
			//	use normal approximation for cube root(chi-squared/dof)
			vd = 1.0 - 2.0 / (9.0 * dof) + 80.0 / (2187.0 * dof * dof * dof);
			vb = sqrt(2.0 / (9.0 * dof) - 104.0 / (2187.0 * dof * dof * dof));
		}	//	end of else for if dof < 25
	}		//	end of else for if dof = 12

}

double rand_u01(unsigned int seed[])
{
//	This code is an exact copy of the C code in L'Ecuyer (Operations Research 1999),
//	with double precision seeds replaced with unsigned integer seeds.  Calculations for 
//	integer seeds are performed with 64 bit integers (long long), and results are
//	stored as 32 bit unsigned integers (unsigned int).  Integer calculations are 
//	faster if using 64 bit configuration (x64) in Release mode (O2 optimization)
//	See L’Ecuyer, P. (1999), “Good Parameters and Implementations for Combined 
//	Multiple Recur-sive Random Number Generators,” Operations Research 47(1):159-164

	double f;
	long long lp1, lp2;

	lp1 = seed[1];
	lp2 = seed[2];
	lp1 = ia12 * lp1 - ia13n * lp2;
	lp1 = lp1 % im1;
	if (lp1 < 0) lp1 += im1;

	seed[2] = seed[1]; seed[1] = seed[0];
	seed[0] = unsigned int(lp1);

	lp1 = seed[3];
	lp2 = seed[5];
	lp2 = ia21 * lp1 - ia23n * lp2;
	lp2 = lp2 % im2;
	if (lp2 < 0) lp2 += im2;

	seed[5] = seed[4]; seed[4] = seed[3];
	seed[3] = unsigned int(lp2);

	if (seed[0] <= seed[3]) f = (((seed[0] - seed[3]) + im1) * norm);
	else f = (double(seed[0] - seed[3]) * norm);

	return f;

}

double sninvdev(unsigned int seed[])
{
	//	This function simulates a standard normal by inverting a U(0,1) using the normal 
	//	inverse function method.  The normal inverse calculation below is accurate to 
	//	1.0 x 10-8.  One can alternatively use an erf_inv function (see Boost library)
	//	Simulation of the U(0,1) is a copy of the code in the function above, 
	//	double rand_u01(unsigned int seed[])

	double p, ans;
	long long lp1, lp2;

	lp1 = seed[1];
	lp2 = seed[2];
	lp1 = ia12 * lp1 - ia13n * lp2;
	lp1 = lp1 % im1;
	if (lp1 < 0) lp1 += im1;

	seed[2] = seed[1]; seed[1] = seed[0];
	seed[0] = unsigned int(lp1);

	lp1 = seed[3];
	lp2 = seed[5];
	lp2 = ia21 * lp1 - ia23n * lp2;
	lp2 = lp2 % im2;
	if (lp2 < 0) lp2 += im2;

	seed[5] = seed[4]; seed[4] = seed[3];
	seed[3] = unsigned int(lp2);

	if (seed[0] <= seed[3]) p = (double((seed[0] - seed[3]) + im1) * norm);
	else p = (double(seed[0] - seed[3]) * norm);

//	ans = sqrt(2.0) * boost::math::erf_inv(2.0 * p - 1.0);

	double q, r;
	if (p <= 0.0) {
		ans = -10.0;
	}
	else {
		if (p >= 1.0) ans = 10.0;

		else {
			if (p < 0.02425) {
				q = sqrt(-2.0 * log(p));
				ans = (((((-0.007784894002430293 * q 
					- 0.3223964580411365) * q 
					- 2.400758277161838) * q 
					- 2.549732539343734) * q 
					+ 4.374664141464968) * q 
					+ 2.938163982698783) /
					((((0.007784695709041462 * q 
						+ 0.3224671290700398) * q 
						+ 2.445134137142996) * q 
						+ 3.754408661907416) * q
						+ 1.0);
			}
			else {
				if (p < 0.97575) {
					q = p - 0.5;
					r = q * q;
					ans = (((((-39.69683028665376 * r 
						+ 220.9460984245205) * r 
						- 275.9285104469687) * r 
						+ 138.3577518672690) * r 
						- 30.66479806614716) * r 
						+ 2.506628277459239) * q /
						(((((-54.47609879822406 * r 
							+ 161.5858368580409) * r 
							- 155.6989798598866) * r 
							+ 66.80131188771972) * r 
							- 13.28068155288572) * r + 1.0);
				}
				else {

					q = sqrt(-2.0 * log(1.0 - p));
					ans = -(((((-0.007784894002430293 * q 
						- 0.3223964580411365) * q 
						- 2.400758277161838) * q 
						- 2.549732539343734) * q 
						+ 4.374664141464968) * q 
						+ 2.938163982698783) /
						((((0.007784695709041462 * q 
							+ 0.3224671290700398) * q
							+ 2.445134137142996) * q 
							+ 3.754408661907416) * q + 1.0);
				}
			}
		}
	}

	return ans;

}

double stdNormalInverse(double p) 
{
	double q, r, ans;
	if (p <= 0.0) {
		ans = -10.0;
	}
	else {
		if (p >= 1.0) ans = 10.0;
		else {
			if (p < 0.02425) {
				q = sqrt(-2.0 * log(p));
				ans = (((((-0.007784894002430293 * q
					- 0.3223964580411365) * q
					- 2.400758277161838) * q
					- 2.549732539343734) * q
					+ 4.374664141464968) * q
					+ 2.938163982698783) /
					((((0.007784695709041462 * q
						+ 0.3224671290700398) * q
						+ 2.445134137142996) * q
						+ 3.754408661907416) * q
						+ 1.0);
			}
			else {
				if (p < 0.97575) {
					q = p - 0.5;
					r = q * q;
					ans = (((((-39.69683028665376 * r
						+ 220.9460984245205) * r
						- 275.9285104469687) * r
						+ 138.3577518672690) * r
						- 30.66479806614716) * r
						+ 2.506628277459239) * q /
						(((((-54.47609879822406 * r
							+ 161.5858368580409) * r
							- 155.6989798598866) * r
							+ 66.80131188771972) * r
							- 13.28068155288572) * r + 1.0);
				}
				else {
					q = sqrt(-2.0 * log(1.0 - p));
					ans = -(((((-0.007784894002430293 * q
						- 0.3223964580411365) * q
						- 2.400758277161838) * q
						- 2.549732539343734) * q
						+ 4.374664141464968) * q
						+ 2.938163982698783) /
						((((0.007784695709041462 * q
							+ 0.3224671290700398) * q
							+ 2.445134137142996) * q
							+ 3.754408661907416) * q + 1.0);
				}
			}
		}
	}
	return ans;
}

double rand_log_u_compl(unsigned int seed[])
{
	//	This function simulates U(0,1) and returns log(1 - U)
	//	U(0,1) simulation uses an exact copy of the C code in L'Ecuyer 

	double u_compl, u, answer;
	long long lp1, lp2;

	lp1 = seed[1];
	lp2 = seed[2];
	lp1 = ia12 * lp1 - ia13n * lp2;
	lp1 = lp1 % im1;
	if (lp1 < 0) lp1 += im1;

	seed[2] = seed[1]; seed[1] = seed[0];
	seed[0] = unsigned int(lp1);

	lp1 = seed[3];
	lp2 = seed[5];
	lp2 = ia21 * lp1 - ia23n * lp2;
	lp2 = lp2 % im2;
	if (lp2 < 0) lp2 += im2;

	seed[5] = seed[4]; seed[4] = seed[3];
	seed[3] = unsigned int(lp2);

	if (seed[0] <= seed[3]) {
		u = (double((seed[0] - seed[3]) + im1) * norm);
		u_compl = u_compl = (double(seed[3] - seed[0]) * norm);
	}
	else {
		u = (double(seed[0] - seed[3]) * norm);
		u_compl = (double(im1 - (seed[0] - seed[3])) * norm);
	}

	if (u < 1.0e-04) answer = -u * (1.0 + 0.5*u + u * u / 3.0 + u * u*u / 4.0 + u * u*u*u / 5.0);
	else {
		if (u < 0.50) answer = log(1 - u);
		else answer = log(u_compl);
	}
	return answer;
}

double chi2_dev(double dof, unsigned int seed[])
{
//	This function simulates a chi squared variate with dof degress of freedom
//	The code below is based on simulation of a gamma variate with shape parameter
//	alpha = dof/2 and beta = 1
//	Simulation method depends on the degrees of freedom parameter, alpha = dof/2

	double term, tem, ans, U1, U2, Bx;
	double b, d, y, p, q, dof1, frtest;
	int i, n;
	const double h_pi = 3.14159265358979;
	const double h_dsqr2 = 0.707106781186547;
	const double h_dsqrpi = 0.564189583547756;

	if (dof <= 0.0) return -9999.99;
	if (dof <= 2.0) {
		if (fabs(dof - 1.0) < 1.0e-06) {
			tem = sninvdev(seed);
			ans = tem * tem;
		}
		else if (fabs(dof - 2.0) < 1.0e-06) ans = -2.0 * rand_log_u_compl(seed);
		else {
		//	Use Johnk's acceptance/rejection method that combines
		//	the Beta distribution with the exponential dostribution
		//	2 simuations with minimum efficiency .81 + 2 final simulations
		//  See Wallace, N.D. (1974), “Computer Generation of Gamma Random Variates
		//	with Non-integral Shape Parameters,” Communications of the Association
		//	for Computing Machinery (ACM) 17 (12), December 1974
			term = 0.5 * dof;
			U1 = rand_u01(seed);
			U2 = rand_u01(seed);
			Bx = pow(U1, 1.0 / term);
			tem = Bx + pow(U2, 1.0 / (1.0 - term));
			while (tem > 1.0) {
				U1 = rand_u01(seed);
				U2 = rand_u01(seed);
				Bx = pow(U1, 1.0 / term);
				tem = Bx + pow(U2, 1.0 / (1.0 - term));
			}
			U1 = rand_u01(seed);
			U2 = rand_u01(seed);
			ans = -log(U1 * U2) * Bx * 2.0;
		}
	}	//	end of if (dof <= 2)
	else {
		if (dof <= 12.0) {
			//	Use Wallace rejection method with 2 exponential simulations
			//	Exponentials for dof = 2, 4, 6, 8, 10, 12
			//	Gamma alpha, a = 1, 2, 3, 4, 5
			//  Up to 6 simulations with efficiency of acceptance/rejection
			//  See Wallace, N.D. (1974), “Computer Generation of Gamma Random Variates
			//	with Non-integral Shape Parameters,” Communications of the Association
			//	for Computing Machinery (ACM) 17 (12), December 1974
			double a = dof / 2;
			bool is_int, is_half_int;
			dof1 = floor(a);
			n = int(dof1);
			is_int = (fabs(dof1 - a) < 1.0e-03);
			is_half_int = (fabs(dof1 - (a + 0.5)) < 1.0e-03);
			if (is_int) {
			//	simulate using the sum of n exponentials
				y = 0.0;
				for (i = 1; i <= n; i++) y += rand_log_u_compl(seed);
				y = -2.0 * y;
				ans = y;
			}
			else if (is_half_int) {
			//	simulate using the sum of n exponentials + chi squared with df=1
				y = 0.0;
				for (i = 1; i <= n; i++) y += rand_log_u_compl(seed);
				y = -2.0 * y;
				tem = sninvdev(seed);
				y += tem * tem;
				ans = y;
			}
			else {
			//	simulate with the probability rejection method of Wallace 1974
				q = a - dof1;
				p = 1.0 - q;
				y = 0.0;
				for (i = 1; i <= n; i++) y -= rand_log_u_compl(seed);
				U1 = rand_u01(seed);
				if (U1 > p) y -= rand_log_u_compl(seed);
				U2 = rand_u01(seed);
				tem = y / dof1;
				frtest = pow(tem, q) / (1.0 + (tem-1.0)*q);
				while (U2 > frtest) {
					y = 0.0;
					for (i = 1; i <= n; i++) y -= rand_log_u_compl(seed);
					U1 = rand_u01(seed);
					if (U1 > p) y -= rand_log_u_compl(seed);
					U2 = rand_u01(seed);
					tem = y / dof1;
					frtest = pow(tem, q) / (1.0 + (tem - 1.0) * q);
				}
				ans = 2.0 * y;
			}
		}	//	end of if (dof <= 12)
		else {
			if (dof <= 25) {
			//	Cheng-Feast (1979) GKM1 ratio of uniform random numbers with acceptance/rejection
			//  See Cheng, R.C.H., and G. M. Feast (1979), “Some Simple Gamma Variate Generators,”
			//  Journal of the Royal Statistical Society, Series C (Applied Statistics), Vol. 28,
			//  No. 3 (1979):  290-295.
			//	consider GMK2 method in Cheng  Feast
				dof1 = 0.5*dof - 1.0;
				b = (0.5*dof - (1.0 / (3.0*dof))) / dof1;
				p = 2.0 / dof1;
				d = p + 2;
				U1 = rand_u01(seed);
				U2 = rand_u01(seed);
				term = b * U2 / U1;
				int ncheck = 0;
				if ((p*U1 - d + term)*term + 1 <= 0.0) ncheck = 1;
				else if (p*log(U1) - log(term) + term - 1.0 <= 0.0) ncheck = 1;
				while (ncheck < 1) {
					U1 = rand_u01(seed);
					U2 = rand_u01(seed);
					term = b * U2 / U1;
					if ((p*U1 - d + term)*term + 1.0 <= 0.0) ncheck = 1;
					else if (p*log(U1) - log(term) + term - 1.0 <= 0.0) ncheck = 1;
				}
				ans = 2.0*dof1*term;
			}		//	end of else for if (dof <= 25)
			else {
			//	use normal approximation for cube root(chi-squared/dof), Hilferty-Wilson method
			//	Wilson, E. B., and M. M. Hilferty. 1931. “The Distribution of Chi-Square.”
			//	Proceedings of the National Academy of Sciences. USA 17 (12): 684–688
				d = 1.0 - 2.0 / (9.0 * dof) + 80.0 / (2187.0 * dof * dof * dof);
				b = sqrt(2.0 / (9.0 * dof) - 104.0 / (2187.0 * dof * dof * dof));
				tem = d + b * sninvdev(seed);
				ans = dof * tem * tem * tem;
				if (ans < 0.0) ans = 0.0;
			}	//	end of else for if (dof <= 25)	
		}	//	end of else for if (dof <= 12)
	}	//	end of else for if (dof <= 2)

	return ans;
}

double chi2dev_set(double dof, unsigned int seed[], double b, double d, double p, 
	int n, bool is_int, bool is_half_int)
{
	//	This function simulates a chi squared variate with dof degress of freedom
	//	The code below is based on simulation of a gamma variate with shape parameter
	//	alpha = dof/2 and beta = 1.  Fixed parameters are set and passed to this function
	//	The fixed parameters are set in void setchi2params(double dof_input, double& vb,
	//	double& vd, double& vp, int& nvexp, bool& is_int, bool& is_half_int)

	double term, tem, ans, U1, U2, Bx;
	double dof2, y, q, frtest;
	int i;
	const double h_pi = 3.14159265358979;
	const double h_dsqr2 = 0.707106781186547;
	const double h_dsqrpi = 0.564189583547756;

	if (dof <= 0.0) return -9999.99;

	if (dof <= 2.0) {
		if (fabs(dof-1.0) < 1.0e-04) {
			tem = sninvdev(seed);
			ans = tem * tem;
		}
		else if (fabs(dof-2.0) < 1.0e-04) ans = -2.0 * rand_log_u_compl(seed);
		else {
			//	Use Johnk's acceptance/rejection method that combines
			//	the Beta distribution with the exponential distribution
			//	2 simuations with minimum efficiency .81 + 2 final simulations
			//  See Wallace, N.D. (1974), “Computer Generation of Gamma Random Variates
			//	with Non-integral Shape Parameters,” Communications of the Association
			//	for Computing Machinery (ACM) 17 (12), December 1974
			term = 0.5 * dof;
			U1 = rand_u01(seed);
			U2 = rand_u01(seed);
			Bx = pow(U1, 1.0 / term);
			tem = Bx + pow(U2, 1.0 / (1.0 - term));
			while (tem > 1.0) {
				U1 = rand_u01(seed);
				U2 = rand_u01(seed);
				Bx = pow(U1, 1.0 / term);
				tem = Bx + pow(U2, 1.0 / (1.0 - term));
			}
			U1 = rand_u01(seed);
			U2 = rand_u01(seed);
			ans = -log(U1 * U2) * Bx * 2.0;
		}
	}	//	end of if (dof <= 2)

	else {
		if (dof <= 12.0) {
			//	Use Wallace rejection method with 2 exponential simulations
			//	Exponentials for dof = 2, 4, 6, 8, ...
			//	Gamma alpha, a = 1, 2, 3, 4, 5
			//  Up to 6 simulations with efficiency of acceptance/rejection
			//  See Wallace, N.D. (1974), “Computer Generation of Gamma Random Variates
			//	with Non-integral Shape Parameters,” Communications of the Association
			//	for Computing Machinery (ACM) 17 (12), December 1974
			if (is_int) {
			//	simulate using the sum of n exponentials
				y = 0.0;
				for (i = 1; i <= n; i++) y += rand_log_u_compl(seed);
				y = -2.0 * y;
				ans = y;
			}
			else if (is_half_int) {
			//	simulate using the sum of n exponentials + chi squared with df=1
				y = 0.0;
				for (i = 1; i <= n; i++) y += rand_log_u_compl(seed);
				y = -2.0 * y;
				tem = sninvdev(seed);
				y += tem * tem;

				ans = y;
			}
			else {
			//	simulate with the probability rejection method of Wallace
				dof2 = floor(dof / 2) ;
				q = 1.0 - p;
				y = 0.0;
				for (i = 1; i <= n; i++) y -= rand_log_u_compl(seed);
				U1 = rand_u01(seed);
				if (U1 > p) y -= rand_log_u_compl(seed);
				U2 = rand_u01(seed);
				tem = y / dof2;
				frtest = pow(tem, q) / (1.0 + (tem - 1.0) * q);
				while (U2 > frtest) {
					y = 0.0;
					for (i = 1; i <= n; i++) y -= rand_log_u_compl(seed);
					U1 = rand_u01(seed);
					if (U1 > p) y -= rand_log_u_compl(seed);
					U2 = rand_u01(seed);
					tem = y / dof2;
					frtest = pow(tem, q) / (1.0 + (tem - 1.0) * q);
				}
				ans = 2.0 * y;
			}
		}	//	end of if (dof <= 12)
		else {
			if (dof <= 25) {
			//	Cheng-Feast (1979) GKM1 ratio of uniform random numbers with acceptance/rejection
			//  See Cheng, R.C.H., and G. M. Feast (1979), “Some Simple Gamma Variate Generators,”
			//  Journal of the Royal Statistical Society, Series C (Applied Statistics), Vol. 28,
			//  No. 3 (1979):  290-295.
			//	consider GMK2 method in Cheng  Feast
				U1 = rand_u01(seed);
				U2 = rand_u01(seed);
				term = b * U2 / U1;
				int ncheck = 0;
				if ((p * U1 - d + term) * term + 1 <= 0.0) ncheck = 1;
				else if (p * log(U1) - log(term) + term - 1.0 <= 0.0) ncheck = 1;
				while (ncheck < 1) {
					U1 = rand_u01(seed);
					U2 = rand_u01(seed);
					term = b * U2 / U1;
					if ((p * U1 - d + term) * term + 1.0 <= 0.0) ncheck = 1;
					else if (p * log(U1) - log(term) + term - 1.0 <= 0.0) ncheck = 1;
				}
				ans = 2.0 * (0.5 * dof - 1.0) * term;
			}		//	end of else for if (dof <= 25)
			else {
			//	use normal approximation for cube root(chi-squared/dof), Hilferty-Wilson method
			//	Wilson, E. B., and M. M. Hilferty. 1931. “The Distribution of Chi-Square.”
			//	Proceedings of the National Academy of Sciences. USA 17 (12): 684–688
				tem = d + b * sninvdev(seed);
				ans = dof * tem * tem * tem;
				if (ans < 0.0) ans = 0.0;
			}	//	end of else for if (dof <= 20)	
		}	//	end of else for if (dof <= 12)
	}	//	end of else for if (dof <= 2)

	return ans;
}

//	need to optimize by starting at largest value for weights
double rn_poisson(double u, double x, double log_Factorial[], double wgt[])
{
	//	x is the Poisson lambda parameter
	//	using MRG32k3a, the smallest value of u is 2.328306549295728e-10
	//	For x > 200, this function is not called, and the normal approximation
	//	for the Poisson distribution is used instead

	int i, j, nparm;
	double term, sum, f, apperr;

	nparm = int(floor(x));
	if (nparm > 50) {
		apperr = 1.0;
		term = exp(nparm * log(x) - x - log_Factorial[nparm]);
		wgt[nparm] = term;
		j = nparm - 1;
		while ((apperr > 1.0e-12) && (j >= 0)) {
			wgt[j] = wgt[j+1] * (double(j) + 1.0) / x;
			apperr = wgt[j];
			j -= 1;
		}
		i = j + 1;
		sum = wgt[i];
		if (sum >= u) f = i;
		else {
			while ((sum < u) && (i < nparm)) {
				i += 1;
				sum += wgt[i];
			}
			if (sum > u) f = i;
			else {
				while ((sum < u) && (i < 400)) {
					i += 1;
					term = term * x / i;
					sum = sum + term;
				}
				f = i;
			}
		}
	}	//	end of if (nparm > 50)
	else {
		term = exp(-x);
		i = 0;
		sum = term;
		if (sum >= u) f = 0.0;
		else {
			while ((sum < u) && (i < 400)) {
				i += 1;
				term = term * x / i;
				sum += term;
			}
			f = i;
		}
	}

	return f;
}
