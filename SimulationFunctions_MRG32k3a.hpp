#pragma once
#include "math.h"

const double norm = 2.328306549295728e-10;
const unsigned int im1 = 4294967087;
const unsigned int im2 = 4294944443;
const unsigned int ia12 = 1403580;
const unsigned int ia13n = 810728;
const unsigned int ia21 = 527612;
const unsigned int ia23n = 1370589;

double rand_u01(unsigned int seed[])
{
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

void roll_seed(unsigned int seed[])
{
//	double f;
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

//	if (seed[0] <= seed[3]) f = (((seed[0] - seed[3]) + im1) * norm);
//	else f = (double(seed[0] - seed[3]) * norm);

//	return f;

}

double sninvdev(unsigned int seed[])
{
	//	This code is an exact copy of the C code in L'Ecuyer (Operations Research 1999)
	//	and simulates a standard normal by inverting a U(0,1) using the normal inverse
	//	The normal inverse calculation below is accurate to 1.0 x 10-8
	//	One can alternatively use an erf_inv function (see Boost library)

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

double stdNormalInverse(double p) {
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
	//	The code below is based on simulatin of a gamma variate with shape parameter
	//	alpha = dof/2 and beta = 1
	double term, tem, ans, U1, U2, Bx;
	double b, d, y, p, dof1, frtest;
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
			//	Compute switch probability
			//  Up to 6 simulations with efficiency ? (see Wallace)
			double a = dof / 2;
			bool is_int, is_half_int;
			dof1 = floor(a);
			n = int(dof1);
		//	is_int = (dof1 == a);
		//	is_half_int = is_int ? false : (fabs(dof1 - a) == 0.5f);
			is_int = (fabs(dof1 - a) < 1.0e-02);
			is_half_int = (fabs(dof1 - (a + 0.5)) < 1.0e-02);
			if (is_int) {
			//	simulate using the sum of n exponentials
			//	y = 1.0;
			//	for (i = 1; i <= n; i++) y *= rand_u_compl(seed);
			//	y = -2.0 * log(y);
				y = 0.0;
				for (i = 1; i <= n; i++) y += rand_log_u_compl(seed);
				y = -2.0 * y;
				ans = y;
			}
			else if (is_half_int) {
			//	simulate using the sum of n exponentials + chi squared with df=1
			//	y = 1.0;
			//	for (i = 1; i <= n; i++) y *= rand_u_compl(seed);
			//	y = -2.0 * log(y);
				y = 0.0;
				for (i = 1; i <= n; i++) y += rand_log_u_compl(seed);
				y = -2.0 * y;
				tem = sninvdev(seed);
				y += tem * tem;
				ans = y;
			}
			else {
				dof1 = 2 * dof1;
				p = 1.0 - 0.5 * (dof - dof1);
				//	simulate with the probability rejection method of Wallace
			//	y = 1.0;
			//	for (i = 1; i <= n; i++) y = y * rand_u_compl(seed);
			//	y = -2.0 * log(y);
				y = 0.0;
				for (i = 1; i <= n; i++) y += rand_log_u_compl(seed);
				y = -2.0 * y;
				U1 = rand_u01(seed);
			//	if (U1 > p) y = y - 2.0 * log(rand_u_compl(seed));
				if (U1 > p) y = y - 2.0 * rand_log_u_compl(seed);
				U2 = rand_u01(seed);
				frtest = pow(y / dof1, 0.5 * (dof - dof1)) / (p + (1.0 - p) * y / dof1);
				while (U2 > frtest) {
				//	y = 1.0;
				//	for (i = 1; i <= n; i++) y = y * rand_u_compl(seed);
				//	y = -2.0 * log(y);
					y = 0.0;
					for (i = 1; i <= n; i++) y += rand_log_u_compl(seed);
					y = -2.0 * y;
					U1 = rand_u01(seed);
				//	if (U1 > p) y = y - 2.0 * log(rand_u_compl(seed));
					if (U1 > p) y = y - 2.0 * rand_log_u_compl(seed);
					U2 = rand_u01(seed);
					frtest = pow(y / dof1, 0.5 * (dof - dof1)) / (p + (1.0 - p) * y / dof1);
				}
				ans = y;
			}
		}	//	end of if (dof <= 12)
		else {
			if (dof <= 25) {
			//	Identify simulation method - probability switch method ?
			//	2 simulations withefficiency ?
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
			}		//	end of else for if (dof <= 20)
			else {
			//	use normal approximation for cube root(chi-squared/dof), Hilferty-Wilson method
				d = 1.0 - 2.0 / (9.0 * dof) + 80.0 / (2187.0 * dof * dof * dof);
				b = sqrt(2.0 / (9.0 * dof) - 104.0 / (2187.0 * dof * dof * dof));
				tem = d + b * sninvdev(seed);
				ans = dof * tem * tem * tem;
				if (ans < 0.0) ans = 0.0;
			}	//	end of else for if (dof <= 20)	
		}	//	end of else for if (dof <= 12)
	}	//	end of else for if (dof <= 2)

	return ans;
}

double chi2dev_set(double dof, unsigned int seed[],
	double b, double d, double dof1, double p, int n)
{
	//	This function simulates a chi squared variate with dof degress of freedom
	//	The code below is based on simulatin of a gamma variate with shape parameter
	//	alpha = dof/2 and beta = 1
	double term, tem, ans, U1, U2, Bx;
	double y, frtest;
	int i;
	const double h_pi = 3.14159265358979;
	const double h_dsqr2 = 0.707106781186547;
	const double h_dsqrpi = 0.564189583547756;

	if (dof <= 0.0) return -9999.99;

	if (dof <= 2.0) {
		if (fabs(dof-1.0) < 1.0e-06) {
			tem = sninvdev(seed);
			ans = tem * tem;
		}
		else if (fabs(dof-2.0) < 1.0e-06) ans = -2.0 * rand_log_u_compl(seed);
		else {
			//	Use Johnk's acceptance/rejection method that combines 
			//	the Beta distribution with the exponential dostribution
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
			//	Compute switch probability
			double a = dof / 2;
			bool is_int, is_half_int;
		//	is_int = (dof1 == a);
		//	is_half_int = is_int ? false : (fabs(dof1 - a) == 0.5f);
			//	is_int = (dof1 == a);
			//	is_half_int = is_int ? false : (fabs(dof1 - a) == 0.5f);
			is_int = (fabs(dof1 - a) < 1.0e-02);
			is_half_int = (fabs(dof1 - (a + 0.5)) < 1.0e-02);
			if (is_int) {
				//	simulate using the sum of n exponentials
			//	y = 1.0;
			//	for (i = 1; i <= n; i++) y *= rand_u_compl(seed);
			//	y = -2.0 * log(y);
				y = 0.0;
				for (i = 1; i <= n; i++) y += rand_log_u_compl(seed);
				y = -2.0 * y;
				ans = y;
			}
			else if (is_half_int) {
				//	simulate using the sum of n exponentials + chi squared with df=1
			//	y = 1.0;
			//	for (i = 1; i <= n; i++) y *= rand_u_compl(seed);
			//	y = -2.0 * log(y);
				y = 0.0;
				for (i = 1; i <= n; i++) y += rand_log_u_compl(seed);
				y = -2.0 * y;
				tem = sninvdev(seed);
				y += tem * tem;

				ans = y;
			}
			else {
				//	simulate with the probability rejection method of Wallace
			//	y = 1.0;
			//	for (i = 1; i <= n; i++) y = y * rand_u_compl(seed);
			//	y = -2.0 * log(y);
				y = 0.0;
				for (i = 1; i <= n; i++) y += rand_log_u_compl(seed);
				y = -2.0 * y;
				U1 = rand_u01(seed);
			//	if (U1 > p) y = y - 2.0 * log(rand_u_compl(seed));
				if (U1 > p) y = y - 2.0 * rand_log_u_compl(seed);
				U2 = rand_u01(seed);
				frtest = pow(y / dof1, 0.5 * (dof - dof1)) / (p + (1.0 - p) * y / dof1);
				while (U2 > frtest) {
				//	y = 1.0;
				//	for (i = 1; i <= n; i++) y = y * rand_u_compl(seed);
				//	y = -2.0 * log(y);
					y = 0.0;
					for (i = 1; i <= n; i++) y += rand_log_u_compl(seed);
					y = -2.0 * y;
					U1 = rand_u01(seed);
				//	if (U1 > p) y = y - 2.0 * log(rand_u_compl(seed));
					if (U1 > p) y = y - 2.0 * rand_log_u_compl(seed);
					U2 = rand_u01(seed);
					frtest = pow(y / dof1, 0.5 * (dof - dof1)) / (p + (1.0 - p) * y / dof1);
				}
				ans = y;
			}
		}	//	end of if (dof <= 12)
		else {
			if (dof <= 25) {
				//	Identify simulation method - probability switch method ?
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
				ans = 2.0 * dof1 * term;
			}		//	end of else for if (dof <= 25)
			else {
				//	use normal approximation for cube root(chi-squared/dof), Hilferty-Wilson method
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
	//	Need to modify to optimize calculations
	//	x is the Poisson lambda parameter
	//	using MRG32k3a, the smallest value of u is 2.328306549295728e-10
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

double nc_chi2(double y, double dof, double c, double vrho,
	unsigned int seed[], double b, double d,
	double dof1, double p, int n, double log_Factorial[], double wksp[])
{
	//	This version of the non-central chi-squared simulator simulates z/c plus
	//	calculation of the derivative with respect to the starting y
	double z, parmnc, temsq, tem, ans;
	//		Simulate a non-central chi-squared for y
	if (dof >= 1.0) {
		parmnc = y * c * vrho;
		tem = sninvdev(seed);
		temsq = sqrt(parmnc);
		z = (tem + temsq)*(tem + temsq);
	//	if (dof > 1.0) z = z + chi2dev_set(dof - 1.0, seed, b, d, dof1, p, n);
		if (dof > 1.000001) z = z + chi2_dev(dof - 1.0, seed);
	}
	else {
		parmnc = 0.5*y*c*vrho;
		if (parmnc < 200.0) {
			tem = rand_u01(seed);
			z = rn_poisson(tem, parmnc, log_Factorial, wksp);
		}
		else {
			//	for parmnc/2 >= 200.0, use normal approximation
			tem = sninvdev(seed)*sqrt(parmnc) + parmnc;
			if (tem < 0.0) tem = 0.0;
			z = floor(tem + 0.5);
		}
		tem = dof + 2.0*z;
		if (tem > 0.0) z = chi2_dev(tem, seed);
		else z = 0.0;
	}
	ans = z / c;
	return ans;
}

void setchi2params(double dof_input, double& vb, double& vd, double& vdof1, double& vp, int& nvexp)
{
	//	function to set chi-sqared parameters that will be reused for cases where dof = dof -1
	double dof1, dof;
	int n;

	dof = dof_input - 1.0;
	vb = 0.0;
	vd = 0.0;
	vdof1 = 0.0;
	vp = 0.0;
	nvexp = 0;

	if (dof <= 12.0) {
		//	Use Wallace rejection method with 2 exponential simulations
		//	Exponentials for dof = 2, 4, 6, 8, ...
		//	Gamma alpha, a = 1, 2, 3, 4, 5
		//	Compute switch probability
		double a = dof / 2;
	//	bool is_int, is_half_int;
		dof1 = floor(a);
		n = int(dof1);
	//	is_int = (dof1 == a);
	//	is_half_int = is_int ? false : (fabs(dof1 - a) == 0.5f);
	//	The following line looks wrong for 2 reasons
	//	dof1 = 2.0 * dof1;
		vdof1 = dof;
		nvexp = n;
		vp = 1.0 - 0.5 * (dof - dof1);
	}	//	end of if (dof <= 12)
	else {
		if (dof <= 25) {
			vdof1 = 0.5 * dof - 1.0;
			vb = (0.5 * dof - (1.0 / (3.0 * dof))) / vdof1;
			vp = 2.0 / vdof1;
			vd = vp + 2.0;
			//	need to pass back vdof1 = dof_input - 1.0
			vdof1 = dof;
		}		//		end of else for if (dof <= 20)
		else {
			//	use normal approximation for cube root(chi-squared/dof)
			vd = 1.0 - 2.0 / (9.0 * dof) + 80.0 / (2187.0 * dof * dof * dof);
			vb = sqrt(2.0 / (9.0 * dof) - 104.0 / (2187.0 * dof * dof * dof));
		}	//	end of else for if dof < 25
	}		//	end of else for if dof = 12

}



double chi2_dev_app(double dof, double p_mix, double b_mix, unsigned int seed[])
{
	//	This function simulates a chi squared variate with dof degress of freedom
	//	The code below is based on simulatin of a gamma variate with shape parameter
	//	alpha = dof/2 and beta = 1
	double term, tem, ans, U1, U2;
	double b, d;
//	double y, p, Bx, frtest, dof1;
//	int i, n;
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
			if (dof < 1.0) {
				//	Approximate with mixture of Beta - Chi Squared(1) distributions 
				term = 2.0 / dof;
				U1 = rand_u01(seed);
				if (U1 <= p_mix) {
					tem = U1 / p_mix;
					ans = dof * pow(tem, term);
				}
				else {
					U2 = (U1 - p_mix) / (1.0 - p_mix);
					tem = stdNormalInverse(U2);
					ans = b_mix * tem * tem;
				}

			}
			else {
				U1 = sninvdev(seed);
				U2 = U1 * U1;
				tem = sqrt(dof);
				ans = (dof - tem) + tem * U2;
			}
		}
	}	//	end of if (dof <= 2)
	else {
		if (dof <= 25.0) {
			U1 = rand_u01(seed);
			tem = sqrt(2.0 * dof);
			ans = (dof - tem) - tem * log(U1);
		}
		else {
			//	use normal approximation for cube root(chi-squared/dof), Hilferty-Wilson method
			d = 1.0 - 2.0 / (9.0 * dof) + 80.0 / (2187.0 * dof * dof * dof);
			b = sqrt(2.0 / (9.0 * dof) - 104.0 / (2187.0 * dof * dof * dof));
			tem = d + b * sninvdev(seed);
			ans = dof * tem * tem * tem;
			if (ans < 0.0) ans = 0.0;
		}	//	end of else for if (dof <= 25)
	}	//	end of else for if (dof <= 2)

	return ans;
}

double chi2dev_set_app(double dof, unsigned int seed[], double p_mix, double b_mix, 
	double b, double d)
{
	//	This function simulates a chi squared variate with dof degress of freedom
	//	The code below is based on simulatin of a gamma variate with shape parameter
	//	alpha = dof/2 and beta = 1
	double term, tem, ans, U1, U2;
//	double y, frtest, Bx;
//	int i;
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
			if (dof < 1.0) {
				//	Approximate with mixture of Beta - Chi Squared(1) distributions 
				term = 2.0 / dof;
				U1 = rand_u01(seed);
				if (U1 <= p_mix) {
					tem = U1 / p_mix;
					ans = dof * pow(tem, term);
				}
				else {
					U2 = (U1 - p_mix) / (1.0 - p_mix);
					tem = stdNormalInverse(U2);
					ans = b_mix * tem * tem;
				}
			}
			else {
				U1 = sninvdev(seed);
				U2 = U1 * U1;
				tem = sqrt(dof);
				ans = (dof - tem) + tem * U2;
			}
		}
	}	//	end of if (dof <= 2)
	else {
		if (dof <= 25.0) {
			U1 = rand_u01(seed);
			tem = sqrt(2.0 * dof);
			ans = (dof - tem) - tem * log(U1);
		}
		else {
			//	use normal approximation for cube root(chi-squared/dof), Hilferty-Wilson method
			tem = d + b * sninvdev(seed);
			ans = dof * tem * tem * tem;
			if (ans < 0.0) ans = 0.0;
		}	//	end of else for if (dof <= 25)
	}	//	end of else for if (dof <= 2)

	return ans;
}

double nc_chi2_app(double y, double dof, double c, double vrho, unsigned int seed[], 
	double b, double d, double p_mix, double b_mix, double log_Factorial[], double wksp[])
{
	//	This version of the non-central chi-squared simulator simulates z/c plus
	//	calculation of the derivative with respect to the starting y
	double z, parmnc, temsq, tem, ans;
	//		Simulate a non-central chi-squared for y
	if (dof >= 1.0) {
		parmnc = y * c * vrho;
		tem = sninvdev(seed);
		temsq = sqrt(parmnc);
		z = (tem + temsq) * (tem + temsq);
		if (dof > 1.000001) z = z + chi2dev_set_app(dof - 1.0, seed, p_mix, b_mix, b, d);
	}
	else {
		parmnc = 0.5 * y * c * vrho;
		if (parmnc < 200.0) {
			tem = rand_u01(seed);
			z = rn_poisson(tem, parmnc, log_Factorial, wksp);
		}
		else {
			tem = sninvdev(seed) * sqrt(parmnc) + parmnc;
			if (tem < 0.0) tem = 0.0;
			z = floor(tem + 0.5);
		}
		tem = dof + 2.0 * z;
		if (tem > 0.0) z = chi2_dev_app(tem, p_mix, b_mix, seed);
		else z = 0.0;
	}
	ans = z / c;
	return ans;
}
