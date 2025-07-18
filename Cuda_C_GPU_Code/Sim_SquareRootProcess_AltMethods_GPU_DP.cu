// SimSquareRootProcess_GPU_float.cu : This file contains the 'main' function. Program execution begins and ends there.
//
//	must add boost folder, go to View Propery Pages -> VC++ Directories -> Include Directories
//
//	This code runs both double precision and single precision versions of the functions on GPU
//
//	Version 1
//	Revised November 2024
//

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include "simFunctions_CPU.hpp"
#include "ncChiSquaredSimulation_MRG32k3a_GPU_DP.cuh"
#include <sys/timeb.h>
#include <time.h>
#include <algorithm>    // std::sort
#include <vector>       // std::vector

using namespace std;

int main()
{
	int i, j, k, nSim, nT, nTimeStepsPerDay, nDaysPerTimeUnit;
	double df, tem, tem_df;
	double time0, time1, time2, time3, time4, time5;
	double x0, kappa, theta, sigma, lambda, c, pmean, pvariance;
	unsigned int iseed[6];
	double dt, b, term1_dt, df_term, sig_dt2, temx, dfterm2;
	double adj1, adj2;
	double vb, vd, vp;
	int nvexp;
	bool is_int, is_half_int;
	double mu_Beta, m2_Beta, Beta_power = 0.0, QE_sig_dt2;
	double sim_mean[6], sim_std[6], sim_m4[6], KS, AD, CvM;
	char HeaderLine[500];
	struct _timeb timebuffer;
	_ftime_s(&timebuffer);
	errno_t errcheck;

	int dcount, nDevices, driverVersion = 0, runtimeVersion = 0;
	cudaDeviceProp deviceProp;
	cudaGetDeviceCount(&nDevices);
	printf(" cudaGetDeviceCount = %i \n", nDevices);

	if (!(nDevices > 0 && nDevices <50)) {
		printf(" This computer does not have a GPU for running Cuda code \n");
		printf(" Enter a number to exit \n");
		std::cin >> i;
		exit(2);
	}

	for (i = 0; i < nDevices; i++) {
		cudaGetDeviceProperties(&deviceProp, i);
		printf("\n Device %d: \"%s\"\n", i, deviceProp.name);
		cudaDriverGetVersion(&driverVersion);
		cudaRuntimeGetVersion(&runtimeVersion);
		printf("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n", 
			driverVersion / 1000, (driverVersion % 100) / 10, runtimeVersion / 1000, 
			(runtimeVersion % 100) / 10);
		printf("  CUDA Capability Major/Minor version number:    %d.%d\n", 
			deviceProp.major, deviceProp.minor);
		printf("  Total amount of global memory:                 %.0f MBytes (%llu bytes)\n",
			(float)deviceProp.totalGlobalMem / 1048576.0f, 
			(unsigned long long) deviceProp.totalGlobalMem);
		printf("  GPU Max Clock rate:                            %.0f MHz (%0.2f GHz)\n", 
			deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);

		cudaDeviceGetAttribute(&dcount, cudaDevAttrMultiProcessorCount, i);
		printf("  cudaDeviceAttribute: MultiProcessorCount = %i \n", dcount);
		cudaDeviceGetAttribute(&dcount, cudaDevAttrMaxThreadsPerMultiProcessor, i);
		printf("  cudaDeviceAttribute: MaxThreadsPerMultiProcessor = %i \n", dcount);
	}
	printf(" \n");

//	default values for tests runs
	iseed[0] = 2739458;
	iseed[1] = 837566;
	iseed[2] = 6626355;
	iseed[3] = 4927348;
	iseed[4] = 1939394;
	iseed[5] = 9834488;

	kappa = 0.25;
	theta = 0.04;
	sigma = 0.5;
	lambda = -0.125;
	x0 = 0.04;

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

	printf(" Model parameters: kappa %f theta %f sigma %f lambda %f x0 %f  degrees of freedom %f \n", kappa, theta,
		sigma, lambda, x0, 4.0 * kappa * theta / (sigma * sigma));

	std::cout << " Enter nT (integer), the length of the time series in days \n";
	std::cin >> nT;
	std::cout << " Enter the number of independent simulations, nSim (integer) \n";
	std::cin >> nSim;
	std::cout << " Enter the number of days per time unit, nDaysPerTimeUnit (integer), 252 or 365 for 1 year \n";
	std::cin >> nDaysPerTimeUnit;
	std::cout << " Now enter time steps per day \n";
	std::cin >> nTimeStepsPerDay ;
	std::cout << " Now enter a new value for sigma (negative value accepts default) \n";
	std::cin >> tem;
	if (tem > 0.0) sigma = tem;
	std::cout << " Now enter a new value for x0 (negative value accepts default) \n";
	std::cin >> tem;
	if (tem > 0.0) x0 = tem;

	printf(" nT = %i nSim = %i nDaysPerTimeUnit %i \n", nT, nSim, nDaysPerTimeUnit);
	fprintf(fout, " nT = %i nSim = %i nDaysPerTimeUnit %i \n", nT, nSim, nDaysPerTimeUnit);
	printf(" sigma %f  and degrees of freedom %f \n", sigma, 4.0 * kappa * theta / (sigma * sigma));

	printf(" Model parameters: kappa %f theta %f sigma %f lambda %f x0 %f \n", kappa, theta,
		sigma, lambda, x0);
	for (i = 0; i < 6; i++) printf(" %u ", iseed[i]);
	printf(" \n");
	fprintf(fout, " Model parameters: kappa %f theta %f sigma %f lambda %f x0 %f  degrees of freedom %f \n", kappa,
		theta, sigma, lambda, x0, 4.0 * kappa * theta / (sigma * sigma));
	for (i = 0; i < 6; i++) fprintf(fout, " %u ", iseed[i]);
	fprintf(fout, " \n");

	double* Sims = new double[6 * long long(nSim)];
	double* CumDistFcn = new double[6 * long long(nSim)];
	double* h_Sims  = new double[long long(nSim)];
	float* Simsf = new float[4 * long long(nSim)];
	float* h_Simsf = new float[long long(nSim)];
	std::vector<double> sortSims(nSim);

	string modelNames[8];
	modelNames[0] = "Quadratic Beta - NC1 approximation";
	modelNames[1] = "Quadratic Beta Double Exponential ";
	modelNames[2] = "Andersen QE approximation         ";
	modelNames[3] = "simple Euler approximation        ";
	modelNames[4] = "Exact Non-Central Chi Squared     ";
	modelNames[5] = "Daily Exact NC Chi Squared        ";

	string floatmodelNames[5];
	floatmodelNames[0] = "Quadratic Beta - NC1 approximation";
	floatmodelNames[1] = "Quadratic Beta Double Exponential ";
	floatmodelNames[2] = "Andersen QE approximation         ";
	floatmodelNames[3] = "simple Euler approximation        ";
	

//	Set up seedsd for GPU device
	unsigned int* h_seeds = new unsigned int[long long(6 * nSim)];
	unsigned int An1[3][3], An2[3][3];

	unsigned long long lp1, lp2, seed1[3], seed2[3], sseed1[3], sseed2[3];
	const unsigned int im1 = 4294967087;
	const unsigned int im2 = 4294944443;

//	Cuda C device memory
	double *d_Sims;
	cudaMalloc((void**)&d_Sims, nSim * sizeof(double));
	unsigned int *d_seeds;
	cudaMalloc((void**)&d_seeds, 6 * nSim * sizeof(unsigned int));
	float* d_Simsf;
	cudaMalloc((void**)&d_Simsf, nSim * sizeof(float));

	int nSimsPerPath;
	nSimsPerPath = nT;

	SkipAhead_MRG32k3a(nSimsPerPath, An1, An2);
//	set seeds for the start of each path
	for (i = 0; i < 6; i++) {
		h_seeds[i] = iseed[i];
	}
	for (i = 0; i < 3; i++) {
		seed1[i] =iseed[i];
		seed2[i] = iseed[i + 3];
	}
	for (k = 1; k < nSim; k++) {
		for (i = 0; i < 3; i++) {
			sseed1[i] = 0;
			sseed2[i] = 0;
			for (j = 0; j < 3; j++) {
				sseed1[i] += (An1[i][j] * seed1[j]) % im1;
				sseed2[i] += (An2[i][j] * seed2[j]) % im2;
			}
			lp1 = sseed1[i];
			lp1 = lp1 % im1;
			sseed1[i] = lp1;
			lp2 = sseed2[i];
			lp2 = lp2 % im2;
			sseed2[i] = lp2;
		}
		for (i = 0; i < 3; i++) {
			h_seeds[i + k * 6] = unsigned int(sseed1[i]);
			h_seeds[i + 3 + k * 6] = unsigned int(sseed2[i]);
		}
		for (i = 0; i < 3; i++) {
			seed1[i] = sseed1[i];
			seed2[i] = sseed2[i];
		}

	}	//	end of loop on k

	cudaMalloc((void**)&d_seeds, 6 * nSim * sizeof(unsigned int));
	cudaMemcpy(d_seeds, h_seeds, 6 * nSim * sizeof(unsigned int), cudaMemcpyHostToDevice);

//	Test non-central chi squared simulators
//	Set constant parameters for simulations
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
	df_term = (kappa * theta - (sigma * sigma / 4)) * term1_dt;
	tem = (kappa + lambda) * dt;
	if (fabs(tem) > 1.0e-06) c = 4.0 * (kappa + lambda) / (sigma * sigma * (1.0 - b));
	else c = 4.0 / (sigma * sigma * dt * (1.0 - 0.5 * tem + tem * tem / 6.0 - tem * tem * tem / 24.0 + tem * tem * tem * tem / 120));
//	c = 4.0 / (sigma * sigma * term1_dt);

	sig_dt2 = sqrt(1.0 / c);
	dfterm2 = kappa * theta * term1_dt;

	QE_sig_dt2 = sigma * sigma * term1_dt;

	if (df < 1.0) {
		//	adjust df for mixed beta exponential
		tem_df = df / 2.0;
		mu_Beta = tem_df / (tem_df + 2.0);
		m2_Beta = tem_df / (tem_df + 4.0);
		Beta_power = 2.0 / tem_df;
	}

	//	Precalculate variance adjustment terms
	adj1 = 0.5 * kappa * theta * term1_dt;
	adj2 = (kappa * theta - sigma * sigma / 8.0) * term1_dt;

//	set fixed parameters for non central chi squared simulation
	setchi2params(df, vb, vd, vp, nvexp, is_int, is_half_int);

	int n_max = 15, max_iter;
	double half_df = 0.5 * df, denom[16];

	denom[0] = tgamma(0.5 * df) * pow(2.0, 0.5 * df);
	for (i = 1; i <= n_max; i++) denom[i] = denom[i - 1] * 2.0 * (0.5 * df + i - 1);
	//	Precalculate variance adjustment terms
	adj1 = 0.5 * kappa * theta * term1_dt;
	adj2 = (kappa * theta - sigma * sigma / 8.0) * term1_dt;

	int* h_iterCount = new int[long long(nSim)];
	double* d_denom;
	int* d_iterCount;
	cudaMalloc((void**)&d_denom, 16 * sizeof(double));
	cudaMalloc((void**)&d_iterCount, nSim * sizeof(int));
	float* d_denom_float;
	cudaMalloc((void**)&d_denom_float, 16 * sizeof(float));
	float denom_f[16];
	for (i = 0; i < 16; i++) denom_f[i] = float(denom[i]);

	cudaMemcpy(d_denom, denom, 16 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_denom_float, denom_f, 16 * sizeof(float), cudaMemcpyHostToDevice);

	if (df_term >= 0.0) {
		modelNames[0] = "Fast Quadratic Approximation      ";
		floatmodelNames[0] = "Fast Quadratic Approximation      ";
	}

	printf(" df_term %16.10e  df %12.6e \n", df_term, df);

	printf(" Running quadratic approximation, either FastQ or QBNC1 \n");
	if (df_term <= 0.0 && df_term > -1.0e-15) {
		df_term = 0.0;
		printf(" Adjusting df_term to be exactly 0.0 in double precision, %12.6e \n",
			df_term);
	}

	printf("     Sometimes a warm up run for the GPU is necessary in this program \n");

	_ftime64_s(&timebuffer);
	time0 = timebuffer.time + timebuffer.millitm / 1000.0;

//	if (df_term >= 0.0) {
//		runFastQ <<<(1 + nSim / 128), 128 >>> (nSim, nT, d_seeds, d_Sims, x0, b, df_term,
//			sig_dt2, dfterm2, c, adj1, adj2);
//	}
//	else {
		runQBNC1 <<<(1 + nSim / 128), 128 >>> (nSim, nT, d_seeds, d_Sims, x0, b, df_term,
			sig_dt2, dfterm2, c, mu_Beta, m2_Beta, Beta_power, adj1, adj2, QE_sig_dt2);
//	}

	cudaGetLastError();
	cudaDeviceSynchronize();

	//Read GPU simulations back to host CPU 
	cudaMemcpy(h_Sims, d_Sims, nSim * sizeof(double), cudaMemcpyDeviceToHost);

	_ftime64_s(&timebuffer);
	time1 = timebuffer.time + timebuffer.millitm / 1000.0;
	time1 = time1 - time0;

	printf("  Compute time for warm up run %7.3f seconds \n", time1);

	printf(" Running quadratic approximation, either FastQ or QBNC1 \n");

	_ftime64_s(&timebuffer);
	time0 = timebuffer.time + timebuffer.millitm / 1000.0;

	if (df_term >= 0.0) {
		runFastQ << <(1 + nSim / 128), 128 >> > (nSim, nT, d_seeds, d_Sims, x0, b, df_term,
			sig_dt2, dfterm2, c, adj1, adj2);
	}
	else {
		runQBNC1 << <(1 + nSim / 128), 128 >> > (nSim, nT, d_seeds, d_Sims, x0, b, df_term,
			sig_dt2, dfterm2, c, mu_Beta, m2_Beta, Beta_power, adj1, adj2, QE_sig_dt2);
	}
		
	cudaGetLastError();
	cudaDeviceSynchronize();

	//Read GPU simulations back to host CPU 
	cudaMemcpy(h_Sims, d_Sims, nSim * sizeof(double), cudaMemcpyDeviceToHost);

	_ftime64_s(&timebuffer);
	time1 = timebuffer.time + timebuffer.millitm / 1000.0;
	time1 = time1 - time0;

	for (i = 0; i < nSim; i++) {
		Sims[i] = h_Sims[i];
	}

	
	printf(" Running simulation of Quadratic Beta Double Exponential \n");
	//	printf(" Running simulation of Quadratic NC Chi Squared approximation \n");

	_ftime64_s(&timebuffer);
	time0 = timebuffer.time + timebuffer.millitm / 1000.0;

	runQB2Exp << <(1 + nSim / 128), 128 >> > (nSim, nT, d_seeds, d_Sims, x0, b,
		df_term, sig_dt2, dfterm2, c, denom[0], half_df, adj1, adj2, QE_sig_dt2, d_iterCount);

	cudaGetLastError();
	cudaDeviceSynchronize();

	//Read GPU simulations back to host CPU 
	cudaMemcpy(h_Sims, d_Sims, nSim * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_iterCount, d_iterCount, nSim * sizeof(int), cudaMemcpyDeviceToHost);

	_ftime64_s(&timebuffer);
	time3 = timebuffer.time + timebuffer.millitm / 1000.0;
	time3 = time3 - time0;

	max_iter = 0;
	for (i = 0; i < nSim; i++) {
		Sims[nSim + i] = h_Sims[i];
		if (h_iterCount[i] > max_iter) max_iter = h_iterCount[i];
	}

	printf(" maximum iterations with simulateQB2ExpNC %i \n", max_iter);
	fprintf(fout, " maximum iterations with simulateQB2Exp %i \n", max_iter);

	//	Andersen's quadratic exponential approximation
	printf(" Running Andersen's quadratic exponential approximation (QE) \n");

	_ftime64_s(&timebuffer);
	time0 = timebuffer.time + timebuffer.millitm / 1000.0;

	runQE <<<(1 + nSim / 128), 128 >>> (nSim, nT, d_seeds, d_Sims, x0, b, QE_sig_dt2, dfterm2);

	cudaGetLastError();
	cudaDeviceSynchronize();

	//Read GPU simulations back to host CPU 
	cudaMemcpy(h_Sims, d_Sims, nSim * sizeof(double), cudaMemcpyDeviceToHost);

	_ftime64_s(&timebuffer);
	time4 = timebuffer.time + timebuffer.millitm / 1000.0;
	time4 = time4 - time0;

	for (i = 0; i < nSim; i++) {
		Sims[2 * nSim + i] = h_Sims[i];
	}

	printf(" Running simple Euler approximation \n");

	//	Simulate the non central chi squared distribution
	_ftime64_s(&timebuffer);
	time0 = timebuffer.time + timebuffer.millitm / 1000.0;

	runEuler <<<(1 + nSim / 128), 128 >>> (nSim, nT, d_seeds, d_Sims, x0, b, kappa, 
		theta, sigma, term1_dt, dt, dfterm2);

	cudaGetLastError();
	cudaDeviceSynchronize();

	//Read GPU simulations back to host CPU 
	cudaMemcpy(h_Sims, d_Sims, nSim * sizeof(double), cudaMemcpyDeviceToHost);

	_ftime64_s(&timebuffer);
	time2 = timebuffer.time + timebuffer.millitm / 1000.0;
	time2 = time2 - time0;

	for (i = 0; i < nSim; i++) {
		Sims[3 * nSim + i] = h_Sims[i];
	}

	nSimsPerPath = nT * 8;
	SkipAhead_MRG32k3a(nSimsPerPath, An1, An2);
	//	set seeds for the start of each path
	for (i = 0; i < 6; i++) {
		h_seeds[i] = iseed[i];
	}
	for (i = 0; i < 3; i++) {
		seed1[i] = iseed[i];
		seed2[i] = iseed[i + 3];
	}
	for (k = 1; k < nSim; k++) {
		for (i = 0; i < 3; i++) {
			sseed1[i] = 0;
			sseed2[i] = 0;
			for (j = 0; j < 3; j++) {
				sseed1[i] += (An1[i][j] * seed1[j]) % im1;
				sseed2[i] += (An2[i][j] * seed2[j]) % im2;
			}
			lp1 = sseed1[i];
			lp1 = lp1 % im1;
			sseed1[i] = lp1;
			lp2 = sseed2[i];
			lp2 = lp2 % im2;
			sseed2[i] = lp2;
		}
		for (i = 0; i < 3; i++) {
			h_seeds[i + k * 6] = unsigned int(sseed1[i]);
			h_seeds[i + 3 + k * 6] = unsigned int(sseed2[i]);
		}
		for (i = 0; i < 3; i++) {
			seed1[i] = sseed1[i];
			seed2[i] = sseed2[i];
		}
	}	//	end of loop on k

	cudaMemcpy(d_seeds, h_seeds, 6 * nSim * sizeof(unsigned int), cudaMemcpyHostToDevice);

	printf(" Running simulation of nc chi squared distribution over time steps \n");

	//	Simulate the non central chi squared distribution
	_ftime64_s(&timebuffer);
	time0 = timebuffer.time + timebuffer.millitm / 1000.0;

	runNC_Chi << <(1 + nSim / 128), 128 >> > (nSim, nT, d_seeds, d_Sims, x0, df, c,
		b, vb, vd, vp, nvexp, is_int, is_half_int);

	cudaGetLastError();
	cudaDeviceSynchronize();

	//Read GPU simulations back to host CPU 
	cudaMemcpy(h_Sims, d_Sims, nSim * sizeof(double), cudaMemcpyDeviceToHost);

	_ftime64_s(&timebuffer);
	time5 = timebuffer.time + timebuffer.millitm / 1000.0;
	time5 = time5 - time0;

	for (i = 0; i < nSim; i++) {
		Sims[4 * nSim + i] = h_Sims[i];
	}


	printf("\n Runtimes:  QBNC1 approximation = %7.3f seconds,  simple Euler approximation = %7.3f seconds  \n",
		time1, time2);
	printf(" Runtimes:  Quadratic Beta Double Exponential = %7.3f seconds   Andersen QE approximation = %7.3f seconds \n",
		time3, time4);
	printf(" Runtimes:  Exact NC Chi Squared Distribution over Daily Time Steps = %7.3f seconds \n",
		time5);

	fprintf(fout, "\n Runtimes:  QBNC1 approximation = %7.3f seconds,  simple Euler approximation = %7.3f seconds  \n",
		time1, time2);
	fprintf(fout, " Runtimes: Quadratic Beta Double Exponential = %7.3f seconds   Andersen QE approximation = %7.3f seconds \n",
		time3, time4);
	fprintf(fout, " Runtimes:  Exact NC Chi Squared Distribution over Daily Time Steps = %7.3f seconds \n",
		time5);

//	Calculate statistics


//	recalculate b and c, and n.c. chi squared parameters
	dt = double(nT) / (double(nDaysPerTimeUnit) * double(nTimeStepsPerDay));
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


	pmean = b * x0 + kappa * theta * term1_dt;
	pvariance = sig_dt2 * (b * x0 + 0.5 * kappa * theta * term1_dt);
	for (j = 0; j < 5; j++) {
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

	//	Use cdf(boost::math::non_central_chi_squared(df, tem), c * Sims[2 * nSim + i]) to calculate K-S
	//	statistics and other measures of goodness of fit

	temx = sqrt(pvariance);
	tem = pvariance * pvariance;

	printf("\n Analytic Solutions Mean, Std Dev.  :      Mean %f  Std Deviation %12.6e \n",
		pmean, sqrt(pvariance));
	fprintf(fout, " Analytic Solutions Mean, Std,Dev.  :      Mean %f  Std Deviation %12.6e \n",
		pmean, sqrt(pvariance));

	for (j = 0; j < 5; j++) {
		printf(" %s :      Mean %f  Std Deviation %12.6e   t tests for mean %6.2f variance %6.2f \n",
			modelNames[j].data(), sim_mean[j], sim_std[j], (sim_mean[j] - pmean) * sqrt(nSim / temx),
			(sim_std[j] * sim_std[j] - pvariance) * sqrt(nSim / (sim_m4[4] - tem)));
		fprintf(fout, " %s :      Mean %f  Std Deviation %12.6e   t tests for mean %6.2f variance %6.2f \n",
			modelNames[j].data(), sim_mean[j], sim_std[j], (sim_mean[j] - pmean) * sqrt(nSim / temx),
			(sim_std[j] * sim_std[j] - pvariance) * sqrt(nSim / (sim_m4[4] - tem)));
	}

	printf("\n Rerunning nc chi squared simulation using one time step over entire period \n");
	printf("   dt %f  b %f  c %f  check df %f \n", dt, b, c, df);

	//	Simulate the non central chi squared distribution
	_ftime64_s(&timebuffer);
	time0 = timebuffer.time + timebuffer.millitm / 1000.0;

	runNC_Chi_Long <<<(1 + nSim / 128), 128 >>> (nSim, d_seeds, d_Sims, x0, df, c,
		b, vb, vd, vp, nvexp, is_int, is_half_int);

	cudaGetLastError();
	cudaDeviceSynchronize();

	//Read GPU simulations back to host CPU 
	cudaMemcpy(h_Sims, d_Sims, nSim * sizeof(double), cudaMemcpyDeviceToHost);

	_ftime64_s(&timebuffer);
	time5 = timebuffer.time + timebuffer.millitm / 1000.0;
	time5 = time5 - time0;

	for (i = 0; i < nSim; i++) {
		Sims[5 * nSim + i] = h_Sims[i];
	}

	printf(" Run time for simulating exact non-central chi squared over longer time interval %f  seconds\n", time5);
	fprintf(fout, " Run time for simulating exact non-central chi squared over longer time interval %f  seconds\n", time5);

	//	Calculate test statistics
	tem = x0 * c * b;
	printf(" ReCalculating nc chi squared distribution for %f degrees of freedom and non-centrality param %f \n", df, tem);

	for (j = 0; j < 6; j++) {
		for (i = 0; i < nSim; i++) {
			CumDistFcn[j * nSim + i] = cdf(boost::math::non_central_chi_squared(df, tem), c * Sims[j * nSim + i]);
		}
	}

	printf("\n Calculating goodness of fit statistics \n");

	for (j = 0; j < 6; j++) {
		for (i = 0; i < nSim; i++) sortSims[i] = CumDistFcn[j * nSim + i];
		std::sort(sortSims.begin(), sortSims.begin() + nSim);
	//	Calculate statistics
	//	ChiTest = 0.0;
		KS = 0.0;
		AD = 0.0;
		CvM = 0.0;
		for (i = 0; i < nSim; i++) {
			tem = sortSims[i] - double(i) / nSim;
			temx = double(i+1) / nSim - sortSims[i];
			if (tem > KS) KS = tem;
			if (temx > KS) KS = temx;
			tem = (double(i + 1) - 0.5) / nSim - sortSims[i];
			CvM += tem * tem;
			if (sortSims[i] <= 0.0 || sortSims[(nSim-1-i)] <= 0.0) AD += -1.0e08;
			else AD -= (2.0 * double(i + 1) - 1.0) * (log(sortSims[i]) + log(1.0 - sortSims[(nSim-1-i)]));
		}
		AD = AD / nSim - nSim; 
		CvM += 1.0 / (12.0 * nSim);
		printf(" Method %i %s:  KS = %7.4f %12.6e  AD = %8.4f  CvM = %10.4f \n", j, modelNames[j].data(), KS, KS, AD, CvM);
		fprintf(fout, " Method %i %s:  KS = %7.4f %12.6e  AD = %10.4f  CvM = %10.4f \n", j, modelNames[j].data(), KS, KS, AD, CvM);
	}

	//	write output for graphing
	//	this can be edited to choose the number of observations to write to the output file
	printf("  Now writing results to output file, first 20,000 simulations \n");
	for (i = 0; i < 20000; i++) {
		for (k = 0; k < 5; k++) {
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
	delete[] h_Sims;
	delete[] h_seeds;

	cudaFree(d_Sims);
	cudaFree(d_seeds);

	cudaDeviceReset();

}
