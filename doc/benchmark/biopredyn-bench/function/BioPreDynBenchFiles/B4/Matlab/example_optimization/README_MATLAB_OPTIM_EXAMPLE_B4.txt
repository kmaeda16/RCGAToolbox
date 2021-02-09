
This is the README file for the example optimization of benchmark B4 using MATLAB in the BioPreDyn-Bench collection.

Created:	17/06/2014
Last modified:	17/06/2014
COntact:	gingproc@iim.csic.es

MATLAB is a high-level language and interactive environment for numerical computation, visualization, and programming. The MATLAB implementations provide model dynamics and calculation of the objective function value. Please note that some aspects of the optimization problem have to be set by the user (initial point, parameter bounds, solver, etc).

SOFTWARE REQUIREMENTS: 
All of the problems in the BioPreDyn-Bench collection can be run in Linux 64 bit and Windows.
A Matlab installation is required. In Windows, the Matlab version must be 32 bit.

HOW TO PROCEED:
The main file is b4_PE.m. It defines a parameter estimation problem and solves it with one of these two methods: fmincon or eSS. 
Open it and edit it according to your preferences. You can modify this file to test any optimization method you wish. 

CONTENTS:
b4_PE.m			defines and solves a parameter estimation problem
b4_bounds.mat 		stores the parameter bounds and initial guess for the optimization.
b4_obj.m 		calculates the objective function used in the optimizations for a particular parameter vector. This function returns three outputs: the objective function itself, the constraints, and the residuals.
radau5g_b4.mexa64	mex-file that is integrated to obtain the model output (Linux, 64 bit) 
radau5g_b4.mexw32 	mex-file that is integrated to obtain the model output (Windows, 32 bit) 
README_MATLAB_OPTIM_EXAMPLE_B4.txt	this file
