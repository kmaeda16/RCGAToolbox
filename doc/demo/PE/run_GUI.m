% This script starts the graphical user interface (GUI) to run a real-coded
% genetic algorithm to solve an example parameter estimation problem.
% 
% 1. Execute this script to open the "RCGAToolbox Mission Cotrol PE"
%    window.
% 2. Push the [Launch] button at the bottom on the window.
% 3. While RCGAs are running, output files are created in the current
%    working directry.
% 
% ------------------------ Example Kinetic Model ------------------------
% - INITIAL CONDITION
% X1 = 0
% X2 = 0
% 
% - PARAMETERS
% X0 = 0.1
% k1 = 1
% k2 = 1
% k3 = 1
% K2 = 1
% K3 = 1
% rootCompartment = 1
% 
% - VARIABLES
% X12 = X1 + X2
% 
% - REACTIONS
% v1 = k1 * X0
% v2 = k2 * (X1/rootCompartment) / (K2 + (X1/rootCompartment))
% v3 = k3 * (X2/rootCompartment) / (K3 + (X2/rootCompartment))
% 
% - BALANCE
% X1_dot = v1 - v2;
% X2_dot = v2 - v3;
% -----------------------------------------------------------------------


% In some versions of MATLAB, you get warnings when the GUI opens. This is
% a known bug of MATLAB. Turn off the warning using warning('off', ...)
warning('off', 'MATLAB:subscripting:noSubscriptsSpecified');

RCGAToolbox_MissionControl_PE;
