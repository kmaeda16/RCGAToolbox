
This is the README file for the MATLAB implementation of benchmark B5 in the BioPreDyn-Bench collection.

http://www.iim.csic.es/~gingproc/biopredynbench/

Created:	17/06/2014
Last modified:	19/06/2014
COntact:	gingproc@iim.csic.es

MATLAB is a high-level language and interactive environment for numerical computation, visualization, and programming. The MATLAB implementations provide model dynamics and calculation of the objective function value. Please note that some aspects of the optimization problem have to be set by the user (initial point, parameter bounds, solver, etc).

SOFTWARE REQUIREMENTS: 
Matlab 64-Bit under Linux, or Matlab 32-Bit under Windows.
For the example below, you will need the MEIGO toolbox:
http://www.iim.csic.es/~gingproc/meigo.html

HOW TO PROCEED:
The Matlab implementation includes a number of files, which allow the user to set up a customized parameter estimation with any method implemented in Matlab.
To test new optimization methods, select b5_obj.m as their objective function. An example of this is provided in the folder "example_optimization".
The type of cost function used can be changed by modifying the script b5_obj.m, which returns the cost as a variable named 'objective', calculated at the end of the script. 
By changing the weights used to calculate the 'objective', different cost functions can be used (note that in this way it is possible to use the same functions as e.g. COPASI).

EXPECTED RESULTS:
We report here, as a reference, an example of the results one can typically expect to obtain:
Using the eSS method with the default settings, we reached a cost function value (expected value to reach) of Jf = 3.0725*10^3 after approximately 8.8*10^4 evaluations (16 hours of CPU time in a computer with Intel Xeon Quadcore processor, 2.50 GHz).
The objective function Jf is of the log-likelihood type, i.e. the squared difference between data and model outputs, divided by the variance. 
Note that these results may vary from one optimization run to the next, due to the stochastic nature of the algorithm.

CONTENTS:
b5_test.m		carries out a test integration.
b5_obj.m		calculates the objective function used in the optimizations for a particular parameter vector. This function returns three outputs: the objective function itself, the constraints, and the residuals.
b5_data.mat		pseudo-experimental data.
cvodesg_b5.mexa64	mex-file that integrates model with CVODES (Linux, Matlab 64 bit) 
cvodesg_b5.mexw32	mex-file that integrates model with CVODES (Windows, Matlab 32 bit) 
README_MATLAB_B5.txt	this file
