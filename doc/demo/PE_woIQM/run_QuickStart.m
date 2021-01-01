% This script demonstrates how to run a real-coded genetic algorithm to
% estimate model parameters in an example kinetic model. This demonstration
% does not rely on IQM Tools.
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


clear mex;
clear all;
close all;

% ======== Experimental Data ========= %
global experimentaldata;
experimentaldata = [ ... % time X1 X2
    0.000E+00 0.000E+00 0.000E+00;
    1.000E+00 6.449E-02 2.559E-02;
    2.000E+00 9.089E-02 5.838E-02;
    3.000E+00 1.023E-01 8.115E-02;
    4.000E+00 1.072E-01 9.489E-02;
    5.000E+00 1.094E-01 1.026E-01;
    6.000E+00 1.104E-01 1.068E-01;
    7.000E+00 1.108E-01 1.089E-01;
    8.000E+00 1.110E-01 1.100E-01;
    9.000E+00 1.111E-01 1.106E-01;
    1.000E+01 1.111E-01 1.109E-01;
    ];

% ========= Problem Settings ========= %
problem.n_gene = 7; % # Decision Variables
problem.n_constraint = 0; % # Constraints
problem.fitnessfun = @Fitness_Example; % Fitness Function
problem.decodingfun = @Decoding_Example; % Decoding Function

% ========= Option Settings ========== %
opts.interimreportfun = @Interimreportfun_Example; % Interim Report Function

% ========== Executing RCGA ========== %
% Results = RCGA_UNDXMGG(problem,opts); % UNDX/MGG
Results = RCGA_REXstarJGG(problem,opts); % REXstar/JGG


%% Print best
paramnames = Model_Example_odefun('parameters');

fprintf('\n--- Best parameter set (f = %e) ---\n',Results.Best.f);
for i = 1 : length(paramnames)
    fprintf('%s = %e\n',char(paramnames(i)),Results.Best.x(i));
end

