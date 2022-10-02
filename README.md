# SimulateSquareRootProcess
C++ code for simulating alternative methods for square root processes

The alternative simulation methods are covered in the paper, "Simulation of Square Root Processes Revisited," by me, Louis Scott.  The 
simulations methods include (1) the exact simulation of a noncentral chi squared, (2) Andersen's quadratic-exponential (QE) method, and (3) 
the quadratic-beta (QB) method developed in the paper.  The QB method uses a quadratic function of a normal simulation when the degrees of 
freedom parameter is greater than or equal to 1.  If the degrees of freedom are less than 1, it uses a mixture of the beta distribution and 
a quadratic normal.

This repository contains the C++ code for running the 3 simulations methods and the simple Euler approximation, as well as the statistical 
tests for the simulations.  The file ncChiSquaredSimulation_MRG32k3a.hpp is a header file with the code to run the QB, QE, and exact non-
central chi squared simulation.  The file Test_NC_ChiSquaredSimulators.cpp is C++ code with a main program to run the statistical tests in
Windows.  To run the full set of statistics, one must download and include the Boost library.  The Boost library function for the noncentral
chi squared function is used to calculate the cumulative distribution function.  The .sln and .vcxproj are the files that can be used on 
Windows to build the executable in Visual Studio.  The file Test_NC_ChiSquareSimulators_Linux.cpp is the version of the C++ code with a main 
program to run on Linux systems.  The file ncChiSquaredSimulation_MRG32k3a_Linux.hpp is the version of the header file that runs on Linux.
To build on Linux, download the two Linux files (.cpp and .hpp) as well as the .csv parameter file and build.  I use Eclipse.

