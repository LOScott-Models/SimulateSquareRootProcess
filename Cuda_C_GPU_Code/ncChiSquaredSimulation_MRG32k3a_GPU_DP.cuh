//	GPU device functions for simulation of square root processes
//
//	Revised November 2024

#include <stdio.h> 

__device__ inline double rand_u01(unsigned int seed[]);
__device__ inline double sninvdev(unsigned int seed[]);
__device__ inline double rand_log_u_compl(unsigned int seed[]);
__device__ inline double chi2_dev(double dof, unsigned int seed[]);
__device__ inline double chi2dev_set(double dof, unsigned int seed[], double b, double d, double p,
	int n, bool is_int, bool is_half_int);
__device__ inline double rn_poisson(double u, double x);

__device__ double simulateFastQ(double x0, unsigned int seed[], double b, double df_term,
	double sig_dt2, double dfterm2, double c, double adj1, double adj2)
{
	double x, sqrt_term, y, tem, tem_mean;
	//	df_term must be >= 0.0
	tem_mean = b * x0;
	sqrt_term = tem_mean + df_term;
	tem = sig_dt2 * sninvdev(seed);
	y = sqrt(sqrt_term);
	x = (y + tem) * (y + tem);
	//	adjust for small variance bias
	tem = tem_mean + dfterm2;
	x = tem + (x - tem) * sqrt((tem_mean + adj1) / (tem_mean + adj2));
	return x;
}

__device__ double simulateQBNC1(double x0, unsigned int seed[], double b, double df_term,
	double sig_dt2, double dfterm2, double c, double mu_Beta, double m2_Beta,
	double Beta_power, double adj1, double adj2, double QE_sig_dt2)
{
	double x, sqrt_term, y, v, tem, quad_a, quad_b, quad_c, p_mix,
		b_mix, lambda_nc, lambda_cf, tem_mean;

	tem_mean = b * x0;
	sqrt_term = tem_mean + df_term;
	if (df_term >= 0.0) {
		tem = sig_dt2 * sninvdev(seed);
		y = sqrt(sqrt_term);
		x = (y + tem) * (y + tem);
		//	adjust for small variance bias
		tem = tem_mean + dfterm2;
		x = tem + (x - tem) * sqrt((tem_mean + adj1) / (tem_mean + adj2));
	}
	else {
		//	This one works better with Andersen's Quadratic here with psi_c = 2.0
		double QE_m, QE_s2, QE_psi, QE_a, QE_b2, QE_b;
		QE_m = tem_mean + dfterm2;
		QE_s2 = QE_sig_dt2 * (b * x0 + 0.5 * dfterm2);
		QE_psi = QE_s2 / (QE_m * QE_m);
		if (QE_psi <= 2.0) {
			tem = 2.0 / QE_psi;
			QE_b2 = tem - 1.0 + sqrt(tem) * sqrt(tem - 1.0);
			QE_b = sqrt(QE_b2);
			QE_a = QE_m / (1.0 + QE_b2);
			tem = QE_b + sninvdev(seed);
			x = QE_a * tem * tem;
		}
	//	double u = rand_u01(seed);
	//	double tem2 = -1.0;
	//	//	sqrt_term >= 0 is equivalent to lambda_nc >= 1 - dof
	//	if (sqrt_term >= 0.0) {
	//		const double sqrt_2 = 1.4142135623731;
	//		//	double tem2;
	//		//	tem = sig_dt2 * sninvdev(seed);
	//		tem = sig_dt2 * sqrt_2 * boost::math::erf_inv(2.0 * u - 1.0);
	//		y = sqrt(sqrt_term);
	//		x = (y + tem) * (y + tem);
	//		//	adjust for small variance bias where you can!!!!
	//		tem = tem_mean + dfterm2;
	//		tem2 = tem + (x - tem) * sqrt((tem_mean + adj1) / (tem_mean + adj2));
	//		//	if (tem2 > 0.0) x = tem2;
	//	}
	//	if (tem2 > 0.0) {
	//		x = tem2;
	//	}
		else {
			lambda_nc = c * tem_mean;
			//	Simulate a mixture of beta and non-central chi square(1) with moment matching
			//	mean
			y = c * (tem_mean + dfterm2);
			//	variance
			v = c * 4.0 * (tem_mean + 0.5 * dfterm2);
			lambda_cf = (lambda_nc * lambda_nc + 6.0 * lambda_nc + 3.0) /
				((lambda_nc + 1.0) * (lambda_nc + 1.0));
			quad_b = (v + y * y) + m2_Beta - 2.0 * mu_Beta * y * lambda_cf;
			quad_a = mu_Beta * mu_Beta * lambda_cf - m2_Beta;
			quad_c = y * y * (lambda_cf - 1.0) - v;
			double tem2 = 4.0 * quad_a * quad_c / (quad_b * quad_b);
			if (fabs(tem2) > 1.0e-04) {
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
			}
			else {
				//	tem = 0.5 * tem2 - tem2 * tem2 / 8.0 + 3.0 * tem2 * tem2 * tem2 / 48.0 
				//										- 5.0 * tem2 * tem2 * tem2 * tem2 / 128.0;
				tem = 0.5 * tem2 * (1.0 - 0.25 * tem2 * (1.0 - 0.5 * tem2 * (1.0 - 5.0 * tem2 / 8.0)));
				p_mix = -tem * 0.5 / quad_a;
				if (p_mix < 0.0 || p_mix > 1.0) printf(" warning in sqrt(1-tem2):  tem2 %12.6e  quad_a %12.6e  quad_c %12.6e  p_mix %f \n",
					tem2, quad_a, quad_c, p_mix);
				b_mix = (y - p_mix * mu_Beta) / ((1.0 - p_mix) * (lambda_nc + 1.0));
			}
			//	Simulate U(0,1) for mixture
			double u = rand_u01(seed);
			if (u <= p_mix) {
				x = pow(u / p_mix, Beta_power);
			}
			else {
			//	const double sqrt_2 = 1.4142135623731;
			//	tem2 = sqrt(lambda_nc) + sqrt_2 * boost::math::erf_inv(2.0 * (1.0 - u) / (1.0 - p_mix) - 1.0);
				tem2 = sqrt(lambda_nc) + normcdfinv((1.0 - u) / (1.0 - p_mix));
				x = b_mix * tem2 * tem2;
			}
			x = x / c;
		}	//	end of else for if (sqrt_term >= 0 && lambda_nc > 2.0)

	}	//	end of else for if (dfterm >= 0)

	return x;
}

__device__ inline double simulateQE(double x0, unsigned int seed[], double b, double QE_sig_dt2,
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
		QE_b2 = (2.0 / QE_psi) - 1.0 + sqrt(2.0 / QE_psi) * sqrt(2.0 / QE_psi - 1.0);
		QE_b = sqrt(QE_b2);
		QE_a = QE_m / (1.0 + QE_b2);
		tem = QE_b + sninvdev(seed);
		x = QE_a * tem * tem;
	}
	else {
		QE_p = (QE_psi - 1.0) / (QE_psi + 1.0);
		QE_beta = 2.0 / (QE_m * (1.0 + QE_psi));
		tem = rand_u01(seed);
		if (tem <= QE_p) x = 0.0;
		else x = -log((1.0 - tem) / (1.0 - QE_p)) / QE_beta;
	}
	return x;
}

__device__ inline double nc_chi2(double x0, double dof, double c, double vrho, unsigned int seed[],
	double b, double d, double p, int n, bool is_int, bool is_half_int)
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
			z = rn_poisson(tem, parmnc);
		}
		else {
			//	for parmnc/2 >= 200.0, use normal approximation
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

__device__ inline double rand_u01(unsigned int seed[])
{
//	The MRG32k3a random number generator is used to simulate a uniform, U(0,1)
//	This code is a copy of the C code in L'Ecuyer (Operations Research 1999), with
//	double precision seeds replaced with unsigned integer seeds.  Calculations for 
//	integer seeds are performed with 64 bit integers (long long), and results are
//	stored as 32 bit unsigned integers (unsigned int).  Integer calculations are 
//	faster if using 64 bit configuration (x64) in Release mode (O2 optimization)
//	See L’Ecuyer, P. (1999), “Good Parameters and Implementations for Combined 
//	Multiple Recursive Random Number Generators,” Operations Research 47(1):159-164

	const double norm = 2.328306549295728e-10;
	const unsigned int im1 = 4294967087;
	const unsigned int im2 = 4294944443;
	const unsigned int ia12 = 1403580;
	const unsigned int ia13n = 810728;
	const unsigned int ia21 = 527612;
	const unsigned int ia23n = 1370589;

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

__device__ inline double sninvdev(unsigned int seed[])
{
//	This function simulates a standard normal by inverting a U(0,1) using the normal 
//	inverse function method.  The Cuda Math library function normcdfinv(x) is used to
//	calculate the inverse of the standard normal distribution function. 
//	Simulation of the U(0,1) here is a copy of the code in the function above, 
//	double rand_u01(unsigned int seed[])

	const double norm = 2.328306549295728e-10;
	const unsigned int im1 = 4294967087;
	const unsigned int im2 = 4294944443;
	const unsigned int ia12 = 1403580;
	const unsigned int ia13n = 810728;
	const unsigned int ia21 = 527612;
	const unsigned int ia23n = 1370589;

	double p, u, ans;
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
	
	if (p <= 0) {
		ans = -10.0;
	}
	else {
		if (p <= 0.99) ans = normcdfinv(p);
		else {
			if (seed[0] <= seed[3]) u = ((1 + seed[3] - seed[0]) * norm);
			else u = ((im1 + 1 - (seed[0] - seed[3])) * norm);
			if (u <= 0) ans = 10.0;
			else ans = -normcdfinv(u);
		}
	}

	return ans;

}

__device__ inline double rand_log_u_compl(unsigned int seed[])
{
//	This function simulates U(0,1) and returns log(1 - U)
//	Uses the code from double rand_u01(unsigned int seed[])

	const double norm = 2.328306549295728e-10;
	const unsigned int im1 = 4294967087;
	const unsigned int im2 = 4294944443;
	const unsigned int ia12 = 1403580;
	const unsigned int ia13n = 810728;
	const unsigned int ia21 = 527612;
	const unsigned int ia23n = 1370589;

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

__device__ inline double chi2_dev(double dof, unsigned int seed[])
{
//	This function simulates a chi squared variate with dof degress of freedom
//	The code below is based on simulation of a gamma variate with shape parameter
//	alpha = dof/2 and beta = 1
//	Simulation method depends on the degrees of freedom parameter, alpha = dof/2

	double term, tem, ans, U1, U2, Bx;
	double b, d, y, p, q, dof1, frtest;
	int i, n;

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

__device__ inline double chi2dev_set(double dof, unsigned int seed[], double b, double d, double p,
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
__device__ inline double rn_poisson(double u, double x)
{
	//	x is the Poisson lambda parameter
	//	using MRG32k3a, the smallest value of u is 2.328306549295728e-10
	int i;
	double term, sum, f;
	
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

	return f;
}

__device__ inline double calculateF_x2(double F_x1, double temexp)
{
	//	double F_x2 = F_x1 + (1.0 - F_x1) * (1.0 - temexp);
	//	return F_x2;

	return F_x1 + (1.0 - F_x1) * (1.0 - temexp);
}

__device__ inline double calculateBeta2(double beta1, double x1, double x2, double adjMean,
	double F_x1, double temexp, double F_x2)
{
	double tem = adjMean - (1.0 - F_x1) * (-(x2 - x1) * temexp + (1.0 / beta1 + x1) * (1.0 - temexp));
	return 1.0 / (tem / (1.0 - F_x2) - x2);
}

__device__ inline double calculateModelMom2(double beta1, double beta2, double x1, double x2,
	double F_x1, double temexp, double F_x2, double beta_pm2)
{
	return beta_pm2 + (1.0 - F_x1) * (-(x2 - x1) * (x2 - x1) * temexp
		+ 2.0 * (1.0 / beta1 + x1) * (-(x2 - x1) * temexp + (1.0 - temexp) / beta1) + x1 * x1 * (1.0 - temexp))
		+ (1.0 - F_x2) * (2.0 * (1.0 / beta2 + x2) / beta2 + x2 * x2);
}

__device__ inline double calculate_dBeta2_dBeta1(double beta1, double beta2, double x1, double x2,
	double F_x1, double temexp, double F_x2, double dF_x2_dBeta1)
{
	return -beta2 * beta2 * ((1.0 - F_x1) * ((1.0 / beta1 + x2) * (x2 - x1) * temexp - (1.0 - temexp) / (beta1 * beta1))
		+ dF_x2_dBeta1 * (1.0 / beta2 + x2)) / (1.0 - F_x2);
}

__device__ inline double calculate_dMom2_dBeta1(double beta1, double beta2, double x1, double x2,
	double F_x1, double temexp, double F_x2, double dF_x2_dBeta1, double dBeta2_dBeta1)
{
	return (1.0 - F_x1) * ((x2 - x1) * (x2 - x1) * (x2 - x1) * temexp
		+ 2.0 * (1.0 / beta1 + x1) * ((x2 - x1) * (x2 - x1) * temexp - (1.0 - temexp)
			/ (beta1 * beta1) + (x2 - x1) * temexp / beta1) - 2.0 * (-(x2 - x1) * temexp
				+ (1.0 - temexp) / beta1) / (beta1 * beta1) + x1 * x1 * (x2 - x1) * temexp)
		- 2.0 * (1.0 - F_x2) * dBeta2_dBeta1 * (2.0 / beta2 + x2) / (beta2 * beta2)
		- dF_x2_dBeta1 * (2.0 * (1.0 / beta2 + x2) / beta2 + x2 * x2);
}

__device__ inline double simulateQB2Exp(double x0, unsigned int seed[], double b, double half_df,
	double df_term, double sig_dt2, double dfterm2, double c, double denom0,
	double adj1, double adj2, double QE_sig_dt2, int& iterCount)
{
	double u, x, sqrt_term, y, tem, tem_mean;

	iterCount = 0;
	tem_mean = b * x0;
	sqrt_term = tem_mean + df_term;
	if (df_term >= 0.0) {
		tem = sig_dt2 * sninvdev(seed);
		y = sqrt(sqrt_term);
		x = (y + tem) * (y + tem);
		//	adjust for small variance bias
		tem = tem_mean + dfterm2;
		x = tem + (x - tem) * sqrt((tem_mean + adj1) / (tem_mean + adj2));
	}
	else {
		int j;
		double lambda_nc = c * tem_mean;
		//	Andersen's Quadratic here with psi = 2.0 or 1.5
		double QE_m, QE_s2, QE_psi, QE_a, QE_b2, QE_b;
		QE_m = tem_mean + dfterm2;
		QE_s2 = QE_sig_dt2 * (tem_mean + 0.5 * dfterm2);
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
		/*
			//	u = rand_u01(seed);
			//	double tem2 = -1.0;
			//	sqrt_term >= 0 is equivalent to lambda_nc >= 1 - dof
				if (sqrt_term >= 0.0) {
					const double sqrt_2 = 1.4142135623731;
					tem = sig_dt2 * sqrt_2 * boost::math::erf_inv(2.0 * u - 1.0);
					y = sqrt(sqrt_term);
					x = (y + tem) * (y + tem);
					//	adjust for small variance bias where you can!!!!
					tem = tem_mean + dfterm2;
					tem2 = tem + (x - tem) * sqrt((tem_mean + adj1) / (tem_mean + adj2));
				}
				if (tem2 > 0.0) {
					x = tem2;
				}
		*/
		else {
			double A_pdf, x_mean, x1, F_x1;
			u = rand_u01(seed);
			A_pdf = exp(-0.5 * lambda_nc) / denom0;
			//	use c * mean = lambda_nc + df here for x_mean
			x_mean = (lambda_nc + 2.0 * half_df);
			x1 = 0.2;
			F_x1 = (A_pdf / half_df) * pow(x1, half_df);
			if (u <= F_x1 || u < 0.01) {
				//	Approximate pdf and cdf using Beta distribution
				x = pow(u * half_df / A_pdf, 1.0 / half_df);
				x = x / c;
			}
			else {
				//	Simulate using inverse CFD method for double exponential
				double x_var, x2, F_x2, beta1, beta2, temexp, err, errThreshold = 0.001;
				double beta_pmean, beta_pm2, modelMom2, x_mean2;
				double dMom2_dBeta1, dBeta2_dBeta1, dF_x2_dBeta1;
				//	adjust x_mean and x1 by dividing by c
				x_mean = x_mean / c;
				x_mean2 = x_mean * x_mean;
				x1 = x1 / c;
				x_var = QE_sig_dt2 * (tem_mean + 0.5 * dfterm2);
				x2 = x_mean + x1;
				beta_pmean = A_pdf * pow(c * x1, half_df + 1.0) / (c * (half_df + 1.0));
				beta_pm2 = A_pdf * pow(c * x1, half_df + 2.0) / (c * c * (half_df + 2.0));
				beta1 = 1.0 / ((x_mean - beta_pmean) / (1.0 - F_x1) - x1);
				temexp = exp(-beta1 * (x2 - x1));
				F_x2 = calculateF_x2(F_x1, temexp);
				beta2 = calculateBeta2(beta1, x1, x2, x_mean - beta_pmean, F_x1, temexp, F_x2);
				modelMom2 = calculateModelMom2(beta1, beta2, x1, x2, F_x1, temexp, F_x2, beta_pm2);
				err = modelMom2 - x_mean2 - x_var;
				j = 0;
				while (j < 25 && (fabs(err) / x_var) > errThreshold) {
					//	Calculate the derivative of modelMom2 with respect to beta1
					dF_x2_dBeta1 = (1.0 - F_x1) * (x2 - x1) * temexp;
					dBeta2_dBeta1 = calculate_dBeta2_dBeta1(beta1, beta2, x1, x2, F_x1, temexp, F_x2, dF_x2_dBeta1);
					dMom2_dBeta1 = calculate_dMom2_dBeta1(beta1, beta2, x1, x2, F_x1, temexp, F_x2, dF_x2_dBeta1, dBeta2_dBeta1);
				//	adjust beta1
				//	double prev_beta1 = beta1;
					beta1 = beta1 - err / dMom2_dBeta1;
					if (beta1 <= 0.0) {
						printf("  Warning:  next iteration for beta1 is negative %12.6e \n", beta1);
						printf(" %i iteration relative error is %12.6e x_init %12.6e x_mean %12.6e  x_var %12.6e  x1, x2 %12.6e %12.6e previous Betas %12.6e  %12.6e dMom2_dBeta1 %12.6e  dF_x2_dBeta1 %12.6e  dBeta2_dBeta1 %12.6e \n",
							j, err / x_var, x0, x_mean, x_var, x1, x2, beta1, beta2, dMom2_dBeta1, dF_x2_dBeta1, dBeta2_dBeta1);
					}
					//	recalculate temexp, F_x2, beta2
					temexp = exp(-beta1 * (x2 - x1));
					F_x2 = calculateF_x2(F_x1, temexp);
					beta2 = calculateBeta2(beta1, x1, x2, x_mean - beta_pmean, F_x1, temexp, F_x2);
					modelMom2 = calculateModelMom2(beta1, beta2, x1, x2, F_x1, temexp, F_x2, beta_pm2);
					err = modelMom2 - x_mean2 - x_var;
					j += 1;
				}	//	end of while (j < 25 && (fabs(err) / x_var) > errThreshold)
				dF_x2_dBeta1 = (1.0 - F_x1) * (x2 - x1) * temexp;
				dBeta2_dBeta1 = calculate_dBeta2_dBeta1(beta1, beta2, x1, x2, F_x1, temexp, F_x2, dF_x2_dBeta1);
				dMom2_dBeta1 = calculate_dMom2_dBeta1(beta1, beta2, x1, x2, F_x1, temexp, F_x2, dF_x2_dBeta1, dBeta2_dBeta1);
				if (j >= 25) {
					//	printf(" Warning:  hit maximum iterations %i and err is %12.6e \n", j, err);
					printf(" max iterations %i u = %12.6e F(x1) %12.6e relative error is %12.6e x_init %12.6e x_mean %12.6e  x_var %12.6e  x1, x2 %12.6e %12.6e  Betas %12.6e  %12.6e dMom2_dBeta1 %12.6e  dF_x2_dBeta1 %12.6e  dBeta2_dBeta1 %12.6e \n",
						j, u, F_x1, err / x_var, x0, x_mean, x_var, x1, x2, beta1, beta2, dMom2_dBeta1, dF_x2_dBeta1, dBeta2_dBeta1);
					printf(" seeds for MRG: %u %u %u %u %u %u \n", seed[0], seed[1], seed[2], seed[3], seed[4], seed[5]);
				}
				//	Use solutions for beta1 and beta 2 to calculate simulations from double exponential
				if (u <= F_x2) {
					x = x1 - log((1.0 - u) / (1.0 - F_x1)) / beta1;
				}
				else {
					x = x2 - log((1.0 - u) / (1.0 - F_x2)) / beta2;
				}

				iterCount = j;
			}

		}	//	end of else for if (sqrt_term >= 0)

	}	//	end of else for if (dfterm >= 0)

	return x;
}

__global__ void runQB2Exp(int nSim, int nT, unsigned int* d_seeds, double* d_Sims,
		double x0, double b, double df_term, double sig_dt2, double dfterm2, double c,
		double denom0, double half_df, double adj1, double adj2, double QE_sig_dt2, 
		int* d_iterCount)
{
	int iSim;
	//Thread index
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;
	iSim = tid;

	if (iSim < nSim) {
		double x;
		unsigned int seed[6];
		int i, iter;
		//	extract seeds for this path from seeds
		for (i = 0; i < 6; i++) seed[i] = d_seeds[i + iSim * 6];
		x = x0;
		for (i = 0; i < nT; i++) {
			x = simulateQB2Exp(x, seed, b, half_df, df_term, sig_dt2, dfterm2, c, denom0,
				adj1, adj2, QE_sig_dt2, iter);
		}
		d_Sims[iSim] = x;
		d_iterCount[iSim] = iter;
	}
}

__global__ void runQBNC1(int nSim, int nT, unsigned int* d_seeds, double *d_Sims, 
	double x0, double b, double df_term, double sig_dt2, double dfterm2, double c,
	double mu_Beta, double m2_Beta, double Beta_power, double adj1, double adj2, 
	double QE_sig_dt2)
{
	int iSim;
	//Thread index
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;
	iSim = tid;

	if (iSim < nSim) {
		double x;
		unsigned int seed[6];
		int i;
		//	extract seeds for this path from seeds
		for (i = 0; i < 6; i++) seed[i] = d_seeds[i + iSim * 6];
		x = x0;
		for (i = 0; i < nT; i++) {
			x = simulateQBNC1(x, seed, b, df_term, sig_dt2, dfterm2, c, mu_Beta, m2_Beta,
				Beta_power, adj1, adj2, QE_sig_dt2);
		}
		d_Sims[iSim] = x;
	}
}

__global__ void runFastQ(int nSim, int nT, unsigned int* d_seeds, double* d_Sims,
	double x0, double b, double df_term, double sig_dt2, double dfterm2, double c,
	double adj1, double adj2)
{
	int iSim;
	//Thread index
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;
	iSim = tid;

	if (iSim < nSim) {
		double x;
		unsigned int seed[6];
		int i;
		//	extract seeds for this path from seeds
		for (i = 0; i < 6; i++) seed[i] = d_seeds[i + iSim * 6];
		x = x0;
		for (i = 0; i < nT; i++) {
			x = simulateFastQ(x, seed, b, df_term, sig_dt2, dfterm2, c, adj1, adj2);
		}
		d_Sims[iSim] = x;
	}
}


__global__ void runQE(int nSim, int nT, unsigned int* d_seeds, double *d_Sims, 
	double x0, double b, double QE_sig_dt2, double dfterm2)
{
	int iSim;
	//Thread index
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;
	iSim = tid;

	if (iSim < nSim) {
		double x;
		unsigned int seed[6];
		int i;
		//	extract seeds for this path from seeds
		for (i = 0; i < 6; i++) seed[i] = d_seeds[i + iSim * 6];
		x = x0;
		for (i = 0; i < nT; i++) {
			x = simulateQE(x, seed, b, QE_sig_dt2, dfterm2);
		}
		d_Sims[iSim] = x;
	}
}

__global__ void runEuler(int nSim, int nT, unsigned int* d_seeds, double* d_Sims, 
	double x0, double b, double kappa, double theta, double sigma, double term1_dt, 
	double dt, double dfterm2)
{
	int iSim;
	//Thread index
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;
	iSim = tid;

	if (iSim < nSim) {
		double x;
		unsigned int seed[6];
		int i;
		//	extract seeds for this path from seeds
		for (i = 0; i < 6; i++) seed[i] = d_seeds[i + iSim * 6];
		x = x0;
		for (i = 0; i < nT; i++) {
			x = b * x + kappa * theta * term1_dt + sigma * sqrt(x * dt) * sninvdev(seed);
			if (x < 0.0) x = 0.0;
		}
		d_Sims[iSim] = x;
	}
}

__global__ void runNC_Chi(int nSim, int nT, unsigned int* d_seeds, double* d_Sims,
	double x0, double df, double c, double vrho, double vb, double vd, double vp, 
	int nvexp, bool is_int, bool is_half_int)
{
	int iSim;
	//Thread index
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;
	iSim = tid;

	if (iSim < nSim) {
		double x;
		unsigned int seed[6];
		int i;
		//	extract seeds for this path from seeds
		for (i = 0; i < 6; i++) seed[i] = d_seeds[i + iSim * 6];
		x = x0;
		for (i = 0; i < nT; i++) {
			x = nc_chi2(x, df, c, vrho, seed, vb, vd, vp, nvexp, is_int, is_half_int);
		}
		d_Sims[iSim] = x;
	}
}

__global__ void runNC_Chi_Long(int nSim, unsigned int* d_seeds, double* d_Sims,
	double x0, double df, double c, double vrho, double vb, double vd, double vp, 
	int nvexp, bool is_int, bool is_half_int)
{
	int iSim;
	//Thread index
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;
	iSim = tid;

	if (iSim < nSim) {
		double x;
		unsigned int seed[6];
		int i;
		//	extract seeds for this path from seeds
		for (i = 0; i < 6; i++) seed[i] = d_seeds[i + iSim * 6];
		x = nc_chi2(x0, df, c, vrho, seed, vb, vd, vp, nvexp, is_int, is_half_int);
		d_Sims[iSim] = x;
	}
}