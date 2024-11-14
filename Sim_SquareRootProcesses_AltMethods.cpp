// Sim_SquareRootProcesses_AltMethods.cpp : This file contains the 'main' function. 
//	Program execution begins and ends there.
//
//	main program uses Boost Library 
//	must add boost folder, go to View Propery Pages -> VC++ Directories -> Include Directories
// 
//	Version 1, to run simulations and test statistics for alternative methods to 
//		simulate square root processes
// 
//	Version 2, includes additional functions for alternative simulation methods
// 
//	Version 3, edited to include only Fast Quadratic (FastQ), QBNC1, QE, and QB2Exp
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include "ncChiSquaredSimulation_MRG32k3a.hpp"
#include <sys/timeb.h>
#include <time.h>
#include <algorithm>    // std::sort
#include <vector>       // std::vector

using namespace std;

int main()
{
	int iT, jT, i, j, k, nSim, nT, nTimeStepsPerDay, nDaysPerTimeUnit;
	double df, tem, tem_df;
	double time0, time1, time2, time3, time4, time5;
	double x, x0, kappa, theta, sigma, lambda, c, pmean, pvariance;
	unsigned int iseed[6], rdseed[6], tdseed[6];
	double dt, b, term1_dt, df_term, sig_dt2, temx, dfterm2;
	double vb, vd, vp;
	int nvexp, nTests;
	bool is_int, is_half_int;
	double mu_Beta = 0.0, m2_Beta = 0.0, Beta_power = 0.0, QE_sig_dt2;
	double sim_mean[8], sim_std[8], sim_m4[8], KS, AD, CvM;
	char HeaderLine[500];
	struct _timeb timebuffer;
	_ftime_s(&timebuffer);
	errno_t errcheck;

	//	default values for tests runs
	iseed[0] = 394995;
	iseed[1] = 45868;
	iseed[2] = 7594993;
	iseed[3] = 1394995;
	iseed[4] = 245868;
	iseed[5] = 79594993;
	kappa = 0.25;
	theta = 0.04;
	sigma = 0.25;
	lambda = -0.125;
	x0 = 0.04;
	nTests = 5;

	nT = 91;
	nSim = 100000;
	nDaysPerTimeUnit = 365;
	nTimeStepsPerDay = 1;

	FILE* fin;
	FILE* fout;
	FILE* fout2;
	FILE* fout3;

	errcheck = fopen_s(&fout2, "SimulatedData.txt", "w");
	errcheck = fopen_s(&fout, "TestNC_ChiSquared_Simulation.txt", "w");
	errcheck = fopen_s(&fout3, "SimulatedDistributionFunctions.txt", "w");

	errcheck = fopen_s(&fin, "SquareRootProcessParams.csv", "r");
	if (errcheck == 0) {
		fgets(HeaderLine, 200, fin);
		fscanf_s(fin, "%lf,%lf,%lf,%lf,%lf", &kappa, &theta, &sigma, &lambda, &x0);
		for (i = 0; i < 6; i++) fscanf_s(fin, ",%u", &iseed[i]);
		fclose(fin);
	}
	else printf(" Cannot find the parameter input file SquareRootProcessParams.csv, using default values \n");

	printf(" Model parameters: kappa %f theta %f sigma %f lambda %f x0 %f  degrees of freedom %f \n", 
		kappa, theta, sigma, lambda, x0, 4.0 * kappa * theta / (sigma * sigma));

	for (i = 0; i < 6; i++) {
		rdseed[i] = iseed[i];
		tdseed[i] = iseed[i];
	}

	std::cout << " Enter nT (integer), the length of the time series in days \n";
	std::cin >> nT;
	std::cout << " Enter the number of independent simulations, nSim (integer) \n";
	std::cin >> nSim;
	std::cout << " Enter the number of days per time unit, nDaysPerTimeUnit (integer), 252 or 365 for 1 year \n";
	std::cin >> nDaysPerTimeUnit;
	std::cout << " Now enter time steps per day \n";
	std::cin >> nTimeStepsPerDay;
	std::cout << " Now enter a new value for sigma (negative value accepts default) \n";
	std::cin >> tem;
	if (tem > 0.0) sigma = tem;

	std::cout << " Now enter a new value for x0 (negative value accepts default) \n";
	std::cin >> tem;
	if (tem > 0.0) x0 = tem;

	std::cout << " Do you want to run the daily simulation of the exact non-central chi squared?  Enter 1 for Yes, 0 for No \n";
	std::cin >> j;
	if (j == 1) {
		nTests = 6;
		std::cout << " Do you want to run inverse CDF method for the daily simulation of the exact non-central chi squared?  Enter 1 for Yes, 0 for No \n";
		std::cin >> i;
		if (i == 1) nTests = 7;
	}

	printf(" nT = %i nSim = %i nDaysPerTimeUnit %i nTimeStepsPerDay %i \n", nT, nSim, nDaysPerTimeUnit, 
		nTimeStepsPerDay);
	fprintf(fout, " nT = %i nSim = %i nDaysPerTimeUnit %i nTimeStepsPerDay %i \n", nT, nSim, nDaysPerTimeUnit, 
		nTimeStepsPerDay);
	printf(" sigma %f  and degrees of freedom %f \n", sigma, 4.0 * kappa * theta / (sigma * sigma));

	printf(" Model parameters: kappa %f theta %f sigma %f lambda %f x0 %f \n", kappa, theta,
		sigma, lambda, x0);
	for (i = 0; i < 6; i++) printf(" %u ", iseed[i]);
	printf(" \n");
	fprintf(fout, " Model parameters: kappa %f theta %f sigma %f lambda %f x0 %f  degrees of freedom %f \n", 
		kappa, theta, sigma, lambda, x0, 4.0 * kappa * theta / (sigma * sigma));
	for (i = 0; i < 6; i++) fprintf(fout, " %u ", iseed[i]);
	fprintf(fout, " \n");

	double* Sims = new double[8 * long long(nSim)];
	double* CumDistFcn = new double[8 * long long(nSim)];
	std::vector<double> sortSims(nSim);
	string modelNames[8];

	modelNames[0] = "Quadratic Beta - NC1 approximation";
	modelNames[1] = "Quadratic Beta Double Exponential ";
	modelNames[2] = "Andersen QE approximation         ";
	modelNames[3] = "simple Euler approximation        ";
	modelNames[4] = "Exact Non-Central Chi Squared     ";
	modelNames[5] = "Daily Exact NC Chi Squared        ";
	modelNames[6] = "Daily Simulation with Inverse CDF ";

	//	Test non-central chi squared simulators
	//	Set constant parameters for simulations upfront
	dt = (1.0 / (double(nDaysPerTimeUnit) * double(nTimeStepsPerDay)));
	b = exp(-(kappa + lambda) * dt);
	printf("   dt %f  and b %f \n", dt, b);
	tem = (kappa + lambda) * dt;
	if (fabs(tem) < 1.0e-04) {
		term1_dt = dt * (1.0 - 0.5 * tem + tem * tem / 6 - tem * tem * tem / 24 + tem * tem * tem * tem / 120);
	}
	else {
		term1_dt = dt * (1 - exp(-tem)) / tem;
	}

	df = 4.0 * kappa * theta / (sigma * sigma);
	df_term = (kappa * theta - (sigma * sigma / 4.0)) * term1_dt;
	tem = (kappa + lambda) * dt;
	if (fabs(tem) > 1.0e-06) c = 4.0 * (kappa + lambda) / (sigma * sigma * (1.0 - b));
	else c = 4.0 / (sigma * sigma * dt * (1.0 - 0.5 * tem + tem * tem 
			/ 6.0 - tem * tem * tem / 24.0 + tem * tem * tem * tem / 120));

	sig_dt2 = sqrt(1.0 / c);
	dfterm2 = kappa * theta * term1_dt;

	//	set fixed parameters for non central chi squared simulation
	setchi2params(df, vb, vd, vp, nvexp, is_int, is_half_int);

	printf(" df_term %16.10e  df %12.6e \n", df_term, df);
	if (df_term <= 0.0 && df_term > -1.0e-15) {
		df_term = 0.0;
		df = 1.0;
		printf(" Adjusting df_term to be exactly 0.0 in double precision, %12.6e and df is set to %f \n",
			df_term, df);
	}

	//	Upfront calculation for Quadratic
	QE_sig_dt2 = sigma * sigma * term1_dt;

	if (df_term < 0.0) {
	//	adjust df for beta distribution, divide df by 2
		tem_df = df / 2.0;
	//	tem_df = df;
		mu_Beta = tem_df / (tem_df + 2.0);
		m2_Beta = tem_df / (tem_df + 4.0);
		Beta_power = 2.0 / tem_df;
	}

	//	Precalculate variance adjustment terms
	double adj1, adj2;
	adj1 = 0.5 * kappa * theta * term1_dt;
	adj2 = (kappa * theta - sigma * sigma / 8.0) * term1_dt;

	int n_max = 20, iterCount, max_iter;
	double half_df, denom[30];
	half_df = 0.5 * df;
	denom[0] = tgamma(0.5 * df) * pow(2.0, 0.5 * df);
	for (i = 1; i <= n_max; i++) denom[i] = denom[i - 1] * 2.0 * (0.5 * df + i - 1);

	if (df_term >= 0.0) {
		printf(" Running fast quadratic approximation (FastQ) for dof >= 1.0 \n");
		modelNames[0] = "Fast Quadratic Approximation      ";
	}
	else printf(" Running quadratic mixed beta approximation, with non-nentral chi squared 1 (QBNC1) \n");
	
	_ftime64_s(&timebuffer);
	time0 = timebuffer.time + timebuffer.millitm / 1000.0;

	if (df_term >= 0.0) {
		for (i = 0; i < nSim; i++) {
			x = x0;
			for (iT = 0; iT < nT; iT++) {
				for (jT = 0; jT < nTimeStepsPerDay; jT++) {
					x = simulateFastQ(x, rdseed, b, df_term, sig_dt2, dfterm2, c, adj1, adj2);
				}
			}
			Sims[i] = x;
		}
	}
	else {
		for (i = 0; i < nSim; i++) {
			x = x0;
			for (iT = 0; iT < nT; iT++) {
				for (jT = 0; jT < nTimeStepsPerDay; jT++) {
					x = simulateQBNC1(x, rdseed, b, df_term, sig_dt2, dfterm2, c, mu_Beta,
						m2_Beta, Beta_power, adj1, adj2, QE_sig_dt2);
				}
			}
			Sims[i] = x;
		}
	}

	_ftime64_s(&timebuffer);
	time1 = timebuffer.time + timebuffer.millitm / 1000.0;
	time1 = time1 - time0;

	//	QB2Exp Method
	//	reset initial seeds so that Quadratic Beta Double Exponential method uses same seeds as QBNC1 method
	for (i = 0; i < 6; i++) rdseed[i] = iseed[i];

	printf(" Runnung simulation of Quadratric with Beta - Double Exponential Approximation \n");

	max_iter = 0;

	_ftime64_s(&timebuffer);
	time0 = timebuffer.time + timebuffer.millitm / 1000.0;

	for (i = 0; i < nSim; i++) {
		x = x0;
		for (iT = 0; iT < nT; iT++) {
			for (jT = 0; jT < nTimeStepsPerDay; jT++) {
				x = simulateQB2Exp(x, rdseed, b, half_df, df_term, sig_dt2, dfterm2, c, denom[0],
					adj1, adj2, QE_sig_dt2, iterCount);
				if (iterCount > max_iter) max_iter = iterCount;
			}
		}
		Sims[nSim + i] = x;
	}

	_ftime64_s(&timebuffer);
	time2 = timebuffer.time + timebuffer.millitm / 1000.0;
	time2 = time2 - time0;

	printf(" maximum iterations with Quadratric Beta - Double Exponential method %i \n", max_iter);
	fprintf(fout, " maximum iterations with Quadratric Beta - Double Exponential method %i \n", max_iter);


	//	Andersen's quadratic exponential (QE) approximation
	//	reset initial seeds so that Andersen's method uses same seeds as QB method
	for (i = 0; i < 6; i++) rdseed[i] = iseed[i];

	printf(" Running Andersen's quadratic exponential approximation (QE) \n");

	//	reset initial seeds so that Andersen's method uses same seeds as QB method
	for (i = 0; i < 6; i++) rdseed[i] = iseed[i];

	_ftime64_s(&timebuffer);
	time0 = timebuffer.time + timebuffer.millitm / 1000.0;

	for (i = 0; i < nSim; i++) {
		x = x0;
		for (iT = 0; iT < nT; iT++) {
			for (jT = 0; jT < nTimeStepsPerDay; jT++) {
				x = simulateQE(x, rdseed, b, QE_sig_dt2, dfterm2);
			}
		}
		Sims[2 * nSim + i] = x;
	}

	_ftime64_s(&timebuffer);
	time3 = timebuffer.time + timebuffer.millitm / 1000.0;
	time3 = time3 - time0;

	//	Euler approximation
	printf(" Running simple Euler approximation \n");
	//	reset initial seeds so that Euler method uses same seeds as QBNC1 method
	for (i = 0; i < 6; i++) rdseed[i] = iseed[i];

	//	Simulate the non central chi squared distribution
	_ftime64_s(&timebuffer);
	time0 = timebuffer.time + timebuffer.millitm / 1000.0;

	for (i = 0; i < nSim; i++) {
		x = x0;
		for (iT = 0; iT < nT; iT++) {
			for (jT = 0; jT < nTimeStepsPerDay; jT++) {
				//	x = (kappa * theta - (kappa + lambda) * x) * dt + sigma * sqrt(x * dt) * sninvdev(rdseed);
				//	instead, use exact discrete time solution for mean
				x = b * x + kappa * theta * term1_dt + sigma * sqrt(x * dt) * sninvdev(rdseed);
				if (x < 0.0) x = 0.0;
			}
		}
		Sims[3 * nSim + i] = x;
	}

	_ftime64_s(&timebuffer);
	time4 = timebuffer.time + timebuffer.millitm / 1000.0;
	time4 = time4 - time0;

	//	The following are used in the simulation of the exact non-central chi squared distribution
	double log_Factorial[200], wksp[200];
	//	Initialize the array for log_Factorial
	log_Factorial[0] = 0.0;
	for (i = 1; i < 200; i++) log_Factorial[i] = log_Factorial[i - 1] + log(double(i));

	printf("\n Runtimes: %s = %7.3f seconds,  %s = %7.3f seconds  \n",
		modelNames[0].data(), time1, modelNames[1].data(), time2);
	printf(" Runtimes: %s = %7.3f seconds   %s = %7.3f seconds      \n",
		modelNames[2].data(), time3, modelNames[3].data(), time4);

	fprintf(fout, "\n Runtimes: %s = %7.3f seconds,  %s = %7.3f seconds  \n",
		modelNames[0].data(), time1, modelNames[1].data(), time2);
	fprintf(fout, " Runtimes: %s = %7.3f seconds   %s = %7.3f seconds      \n",
		modelNames[2].data(), time3, modelNames[3].data(), time4);

	if (nTests > 5) {
		//	set fixed parameters for non central chi squared simulation
		setchi2params(df, vb, vd, vp, nvexp, is_int, is_half_int);

		printf(" Running simulation of nc chi squared distribution over time steps; this can be very slow, please be patient \n");

		//	Simulate the non central chi squared distribution
		_ftime64_s(&timebuffer);
		time0 = timebuffer.time + timebuffer.millitm / 1000.0;

		for (i = 0; i < nSim; i++) {
			x = x0;
			for (iT = 0; iT < nT; iT++) {
				for (jT = 0; jT < nTimeStepsPerDay; jT++) {
					//	Simulate a non-central chi-squared for x using the exact methods	
					x = nc_chi2(x, df, c, b, tdseed, vb, vd, vp, nvexp, is_int,
						is_half_int, log_Factorial, wksp);
				}
			}
			Sims[5 * nSim + i] = x;
		}

		_ftime64_s(&timebuffer);
		time5 = timebuffer.time + timebuffer.millitm / 1000.0;
		time5 = time5 - time0;

		printf(" Runtimes: %s = %7.3f seconds \n", modelNames[5].data(), time5);
		fprintf(fout, " Runtimes: %s = %7.3f seconds \n", modelNames[5].data(), time5);

		if (nTests > 6) {
			printf(" Running simulation of exact NC chi squared with inverse CDF method; this can be very slow, please be patient \n");

			double dof = 4.0 * kappa * theta / (sigma * sigma);

			//	Simulate the non central chi squared distribution using gamma CDF inverse method for chi squared
			//	with 2 simulations per time step
			for (i = 0; i < 6; i++) rdseed[i] = iseed[i];
			max_iter = 0;

			_ftime64_s(&timebuffer);
			time0 = timebuffer.time + timebuffer.millitm / 1000.0;

			for (i = 0; i < nSim; i++) {
				x = x0;
				for (i = 0; i < nSim; i++) {
					x = x0;
					for (iT = 0; iT < nT; iT++) {
						for (jT = 0; jT < nTimeStepsPerDay; jT++) {
							x = nc_chi2invdev(x, dof, c, b, rdseed, log_Factorial, wksp);
						}
					}
					Sims[6 * nSim + i] = x;
				}
			}

			_ftime64_s(&timebuffer);
			double time6 = timebuffer.time + timebuffer.millitm / 1000.0;
			time6 = time6 - time0;

			printf(" Runtimes: %s = %7.3f seconds \n", modelNames[6].data(), time6);
			fprintf(fout, " Runtimes: %s = %7.3f seconds \n", modelNames[6].data(), time6);
		}	//	end of if (nTests > 6)
	}	//	end of if (nTests > 5)

	//	recalculate b and c, and n.c. chi squared parameters for 1 time step over entire period
	dt = double(nT) / (double(nDaysPerTimeUnit));
	b = exp(-(kappa + lambda) * dt);
	tem = (kappa + lambda) * dt;
	if (fabs(tem) > 1.0e-06) c = 4.0 * (kappa + lambda) / (sigma * sigma * (1.0 - b));
	else c = 4.0 / (sigma * sigma * dt * (1.0 - 0.5 * tem + tem * tem / 6.0 - tem * tem * tem / 24.0 + tem * tem * tem * tem / 120));
	if (fabs(tem) < 1.0e-04) {
		term1_dt = dt * (1.0 - 0.5 * tem + tem * tem / 6 - tem * tem * tem / 24 + tem * tem * tem * tem / 120);
	}
	else {
		term1_dt = dt * (1 - exp(-tem)) / tem;
	}
	sig_dt2 = sigma * sigma * term1_dt;

	setchi2params(df, vb, vd, vp, nvexp, is_int, is_half_int);

	printf("\n Rerunning nc chi squared simulation using one time step over entire period \n");
	printf("   dt %f  b %f  c %f  check df %f \n", dt, b, c, df);

	//	Simulate the non central chi squared distribution
	_ftime64_s(&timebuffer);
	time0 = timebuffer.time + timebuffer.millitm / 1000.0;

	for (i = 0; i < nSim; i++) {
		x = x0;
		x = nc_chi2(x, df, c, b, tdseed, vb, vd, vp, nvexp, is_int, is_half_int, log_Factorial, wksp);
		Sims[4 * nSim + i] = x;
	}

	_ftime64_s(&timebuffer);
	time5 = timebuffer.time + timebuffer.millitm / 1000.0;
	time5 = time5 - time0;

	printf(" Run time for simulating exact non-central chi squared over longer time interval %f  seconds\n\n", time5);
	fprintf(fout, " Run time for simulating exact non-central chi squared over longer time interval %f  seconds\n\n", time5);

	//	Calculate statistics

	pmean = b * x0 + kappa * theta * term1_dt;
	pvariance = sig_dt2 * (b * x0 + 0.5 * kappa * theta * term1_dt);

	for (j = 0; j < nTests; j++) {
		sim_mean[j] = 0.0;
		sim_std[j] = 0.0;
		sim_m4[j] = 0.0;
		for (i = 0; i < nSim; i++) {
			sim_mean[j] += Sims[j * nSim + i];
			tem = (Sims[j * nSim + i] - pmean);
			tem = tem * tem;
			sim_std[j] += tem;
			sim_m4[j] += tem * tem;
		}
		sim_mean[j] = sim_mean[j] / double(nSim);
		sim_std[j] = sim_std[j] / double(nSim);
		sim_std[j] = sqrt(sim_std[j]);
		sim_m4[j] = sim_m4[j] / double(nSim);
	}

	temx = sqrt(pvariance);
	tem = pvariance * pvariance;
	//	Use the exact nc chi squared simulation for the 4th moment
	double ncchi_m4 = sim_m4[4];

	printf(" Analytic Solutions Mean, Std Dev. :  Mean %f  Std Deviation %12.6e \n",
		pmean, sqrt(pvariance));
	fprintf(fout, " Analytic Solutions Mean, Std,Dev. :  Mean %f  Std Deviation %12.6e \n",
		pmean, sqrt(pvariance));

	for (i = 0; i < nTests; i++) {
		printf(" %s:  Mean %f  Std Deviation %12.6e   t tests for mean %6.2f variance %6.2f \n",
			modelNames[i].data(), sim_mean[i], sim_std[i], (sim_mean[i] - pmean) * sqrt(nSim / temx),
			(sim_std[i] * sim_std[i] - pvariance) * sqrt(nSim / (ncchi_m4 - tem)));
		fprintf(fout, " %s:  Mean %f  Std Deviation %12.6e   t tests for mean %6.2f variance %6.2f \n",
			modelNames[i].data(), sim_mean[i], sim_std[i], (sim_mean[i] - pmean) * sqrt(nSim / temx),
			(sim_std[i] * sim_std[i] - pvariance) * sqrt(nSim / (ncchi_m4 - tem)));
	}

	//	Use cdf(boost::math::non_central_chi_squared(df, tem), c * Sims[2 * nSim + i]) to calculate K-S
	//	statistics and other measures of goodness of fit

	//	Calculate test statistics
	tem = x0 * c * b;
	printf("  Calculating empirical distributions for nc chi squared for %f degrees of freedom and non-centrality param %f \n", df, tem);

	for (j = 0; j < nTests; j++) {
		for (i = 0; i < nSim; i++) {
			CumDistFcn[j * nSim + i] = cdf(boost::math::non_central_chi_squared(df, tem), c * Sims[j * nSim + i]);
		}
	}

	printf("  Calculating goodness of fit statistics \n");

	for (j = 0; j < nTests; j++) {
		for (i = 0; i < nSim; i++) sortSims[i] = CumDistFcn[j * nSim + i];
		std::sort(sortSims.begin(), sortSims.begin() + nSim);
		//	Calculate statistics
		KS = 0.0;
		AD = 0.0;
		CvM = 0.0;
		for (i = 0; i < nSim; i++) {
			tem = sortSims[i] - double(i) / nSim;
			temx = double(i + 1) / nSim - sortSims[i];
			if (tem > KS) KS = tem;
			if (temx > KS) KS = temx;
			tem = (double(i + 1) - 0.5) / nSim - sortSims[i];
			CvM += tem * tem;
			if (sortSims[i] <= 0.0 || sortSims[(nSim - 1 - i)] <= 0.0) AD += -1.0e08;
			else AD -= (2.0 * double(i + 1) - 1.0) * (log(sortSims[i]) + log(1.0 - sortSims[(nSim - 1 - i)]));
		}
		AD = AD / nSim - nSim;
		CvM += 1.0 / (12.0 * nSim);
		printf(" Method %i %s:  KS = %7.4f %12.6e  AD = %8.4f  CvM = %10.4f \n", j, modelNames[j].data(), KS, KS, AD, CvM);
		fprintf(fout, " Method %i %s:  KS = %7.4f %12.6e  AD = %10.4f  CvM = %10.4f \n", j, modelNames[j].data(), KS, KS, AD, CvM);
	}

	//	write output for graphing
	printf("  Now writing results to output file, first 50,000 simulations \n");
	for (i = 0; i < 50000; i++) {
		for (k = 0; k < nTests; k++) {
			fprintf(fout2, " %12.6e ", Sims[k * nSim + i]);
			fprintf(fout3, " %12.6e ", CumDistFcn[k * nSim + i]);
		}
		fprintf(fout2, " \n");
		fprintf(fout3, " \n");
	}

	fclose(fout);
	fclose(fout2);
	fclose(fout3);

	delete[] Sims;
	delete[] CumDistFcn;

}