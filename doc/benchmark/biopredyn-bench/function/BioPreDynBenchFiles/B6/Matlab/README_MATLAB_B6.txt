
This is the README file for the MATLAB implementation of benchmark B6 in the BioPreDyn-Bench collection.

http://www.iim.csic.es/~gingproc/biopredynbench/

Created:	17/06/2014
Last modified:	19/06/2014
COntact:	gingproc@iim.csic.es

MATLAB is a high-level language and interactive environment for numerical computation, visualization, and programming. The MATLAB implementations provide model dynamics and calculation of the objective function value. Please note that some aspects of the optimization problem have to be set by the user (initial point, parameter bounds, solver, etc).

SOFTWARE REQUIREMENTS: 
Matlab 64-Bit or 32-Bit under Linux

HOW TO PROCEED:
The Matlab implementation includes a number of files, which allow the user to set up a customized parameter estimation with any method implemented in Matlab.
To test new optimization methods, select b6_obj.m as their objective function. An example of this is provided in the folder "example_optimization".

EXPECTED RESULTS:
We report here, as a reference, an example of the results one can typically expect to obtain:
Using the eSS method with the default settings, we reached a cost function value (expected value to reach) of Jf = 1.0833*10^5 after approximately 2*10^6 evaluations (24 hours of CPU time in a computer with Intel Xeon Quadcore processor, 2.50 GHz).
The objective function Jf is of the weighted least squares type, where the weights are inversely related to the level of expression. 
Note that these results may vary from one optimization run to the next, due to the stochastic nature of the algorithm.

CONTENTS:
b6_test.m		carries out a test integration.
b6_obj.m		calculates the objective function used in the optimizations for a particular parameter vector. This function returns three outputs: the objective function itself, the constraints, and the residuals
b6_dyn.m		returns the system outputs that will be used in the calculation of the objective function 
b6.mexa64		mex-file integrating the model (Linux, Matlab 64 bit) 
b6.mexglx		mex-file integrating the model (Linux, Matlab 32 bit) 
dm_hkgn53_wls_5_003	input file for the mex-files, containing model data
callUnfold.m		Matlab wrapper that calls the unfold script with the appropriate command
unfold			returns the model outputs at a particular timepoint
setParameters.m		Matlab wrapper that calls the setParameters script
setParameters		script that changes the value of the model parameters
README_MATLAB_B6.txt	this file