% This script shows how to perform sensitivity analysis using RCGAToolbox. 
% By performing sensitivity analysis, you can check which parameter 
% significantly affects the objective function (f).


clear mex;
clear all;
close all;

% =============== Model =============== %
% model = 'Model_Example_SBML.xml'; % SBML file (IQM Tools required)
% model = IQMmodel('Model_Example_SBML.xml'); % Creating an IQMmodel object (IQM Tools required)
model = @Model_Example_odefun; % ODE file (IQM Tools format)
% model = 'Model_Example_odefun.m'; % ODE file (IQM Tools format)
% model = 'Model_Example_mex.c'; % C source code (IQM Tools required)
% model = 'Model_Example_mex.mexw64'; % MEX model for Windows
% model = 'Model_Example_mex.mexmaci64'; % MEX model file for macOS
% model = 'Model_Example_mex.mexa64'; % MEX model file for Linux

% ========= Reading Parameter Estimation Results ========== %
% load('Report'); % Report.mat has the structure Results
% param = Results.Best.x; % Set the best parameter set obtained in parameter estimation
% param = model('parametervalues'); % Set the default parameter set
param(1) = 0.1; % X0
param(2) = 1;   % k1
param(3) = 1;   % k2
param(4) = 1;   % k3
param(5) = 1;   % K2
param(6) = 1;   % K3
param(7) = 1;   % rootCompartment

% ========= Reading Measurement File ========== %
% measurement = 'Measurement_Example.xls'; % Measurement File (EXCEL format)
measurement = 'Measurement_Example.csv'; % Measurement File (CSV)

% ============ ODE Solver ============= %
fast_flag = 0; % MATLAB ODEXX
% fast_flag = 1; % SundialsTB CVODE (SundialsTB required)
% fast_flag = 2; % IQM Tools CVODE MEX (IQM Tools required)

% ============ Normalization ============= %
% norm_flag = 0; % No normalization for both f and p, i.e., Sensitivity = ( f_ptb - f_std ) / ( param_ptb - param_std )
norm_flag = 1; % Normalization only for p, i.e., Sensitivity = ( f_ptb - f_std ) / ( param_ptb - param_std ) * param_std
% norm_flag = 2; % Normalization only for f, i.e., Sensitivity = ( f_ptb - f_std ) / ( param_ptb - param_std ) / f_std
% norm_flag = 3; % Normalization for both f and p, i.e., Sensitivity = ( f_ptb - f_std ) / ( param_ptb - param_std ) * param_std / f_std 

% ============ Sensitivity Calculation ============= %
S = RCGAsensitivity(model, measurement, param, [], [], [], [], norm_flag);

fprintf('\n');
disp('Sensitivities of f to parameters:');
for i = 1 : length(param)
    fprintf('%s: %e\n',char(S.Properties.RowNames(i)),table2array(S(i,1)));
end
