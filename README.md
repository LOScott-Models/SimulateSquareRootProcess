# SimulateSquareRootProcess
C++ code for simulating alternative methods for square root processes

The alternative simulation methods are covered in the paper, "Simulation of Square Root Processes Revisited,"
by me, Louis Scott.  The simulations methods include (1) the exact simulation of a noncentral chi squared,
(2) Andersen's quadratic-exponential (QE) method, and (3) the quadratic-beta (QB) method developed in the paper.  
The QB method uses a quadratic function of a normal simulation when the degrees of freedom parameter is greater 
than or equal to 1.  If the degrees of freedom are less than 1, it uses a mixture of the beta distribution and
a quadratic normaal.

This repository contains the C++ code for running the 3 simulations methods and the simple Euler approximation, 
as well as the statistical tests for the simulations.  Executables for Windows and Linux and data input files 
are also included.
