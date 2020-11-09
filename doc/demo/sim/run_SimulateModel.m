% This script demonstrates how to run simulation using RCGAToolbox.


clearvars;

% =============== Model =============== %
% model = 'Model_Example_SBML.xml'; % SBML file (IQM Tools required)
% model = IQMmodel('Model_Example_SBML.xml'); % Creating an IQMmodel object (IQM Tools required)
model = @Model_Example_odefun; % ODE file (IQM Tools format)
% model = 'Model_Example_odefun.m'; % ODE file (IQM Tools format)
% model = 'Model_Example_mex.c'; % C source code (IQM Tools required)
% model = 'Model_Example_mex.mexw64'; % MEX model for Windows
% model = 'Model_Example_mex.mexmaci64'; % MEX model file for macOS
% model = 'Model_Example_mex.mexa64'; % MEX model file for Linux

% =============== Time ================ %
tspan = 0 : 0.1 : 10;

% ========= Initial Condition ========= %
y0(1) = 0; % X1
y0(2) = 0; % X2

% ========= Parameter Values ========== %
param(1) = 0.1; % X0
param(2) = 1;   % k1
param(3) = 1;   % k2
param(4) = 1;   % k3
param(5) = 1;   % K2
param(6) = 1;   % K3
param(7) = 1;   % rootCompartment

% ============ ODE Solver ============= %
fast_flag = 0; % # fast_flag (0: MATLAB ODEXX)
% fast_flag = 1; % # fast_flag (1: SundialsTB CVODE) (SundialsTB required)
% fast_flag = 2; % # fast_flag (2: IQM Tools CVODE MEX) (IQM Tools required)

% ============ Simulation ============= %
[ T, Y ] = RCGAsimulate(model, tspan, y0, param, fast_flag);

% =============== Figure ============== %
figure;
plot(T,Y,'-','LineWidth',2);
legend('X_1','X_2','Location','best');
xlabel('Time');
ylabel('Concentration');
