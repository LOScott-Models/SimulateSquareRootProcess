# SimulateSquareRootProcess
C++ code for simulating alternative methods for square root processes

The alternative simulation methods are covered in the paper, "Alternative Methods for Simulating Square Root Processes."  
The simulations methods include (1) the exact simulation of a noncentral chi squared, (2) Andersen's quadratic-exponential (QE) method, (3) 
the quadratic beta double exponential method (QB2Exp) method, and (4) the quadratic beta non-central chi squared 1 method (QBNC1). The QB2Exp
and QBNC1 methods have been developed in the paper. The QB2Exp and the QBNC1 methods use a quadratic function of a normal simulation when 
the degrees of freedom parameter is greater than or equal to 1.  If the degrees of freedom are less than 1, the QB2Exp method alternates 
between the quadratic method and an approximation of the cumulative distribution function (CDF) using a beta distribution and two exponential
distributions. The QBNC1 method alternates between the quadratic method and a mixture of the beta distribution and the non-central chi squared
distribution with one degree of freedom.

This repository contains the C++ code for running the 4 simulations methods and the simple Euler approximation, as well as the statistical 
tests for the simulations. The file ncChiSquaredSimulation_MRG32k3a.hpp is a header file with the code to run the QB2Exp, QBNC1, QE, and 
exact non-central chi squared simulation. The file Test_NC_ChiSquaredSimulators.cpp is C++ code with a main program to run the statistical 
tests in Windows. To run the full set of statistics, one must download and include the Boost library. The Boost library function for the 
noncentral chi squared function is used to calculate the CDF. The .sln and .vcxproj are the files that can be used on Windows to build 
the executable in Visual Studio.  The file Test_NC_ChiSquareSimulators_Linux.cpp is the version of the C++ code with a main program to run 
on Linux systems.  To build on Linux, download Test_NC_ChiSquareSimulators_Linux.cpp, ncChiSquaredSimulation_MRG32k3a.hpp, and the .csv 
parameter file, and build. I use Eclipse to build the executables on Linux.

