% This script demonstrates how to run simulations using RCGAToolbox.


clearvars;

% =============== Model =============== %
% model = 'ExampleModel.xml'; % SBML File
% model = IQMmodel('ExampleModel.xml'); % Creating IQMmodel
model = 'ExampleModel_odefun.m'; % MATLAB ODE Function File
% model = 'ExampleModel_mex.c'; % C ODE File
% model = 'ExampleModel_mex.mexw64'; % MEX ODE File for Windows
% model = 'ExampleModel_mex.mexmaci64'; % MEX ODE File for macOS
% model = 'ExampleModel_mex.mexa64'; % MEX ODE File for Linux

% =============== Time ================ %
tspan = 0 : 0.2 : 20;

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
% fast_flag = 1; % # fast_flag (1: SundialsTB CVODE)
% fast_flag = 2; % # fast_flag (2: IQMTools CVODE MEX)

% ============ Simulation ============= %
[ T, Y ] = RCGAsimulate(model, tspan, y0, param, fast_flag);

% =============== Figure ============== %
figure;
plot(T,Y,'-','LineWidth',2);
legend('X_1','X_2','Location','best');
xlabel('Time');
ylabel('Concentration');
