% This script demonstrates how to run simulations using RCGAToolbox.


clearvars;

% =============== Model =============== %
% modelfile = 'modelExample.xml'; % SBML File
% modelfile = IQMmodel('modelExample.xml'); % Creating IQMmodel
modelfile = 'modelExample_odefun.m'; % MATLAB ODE Function File
% modelfile = 'modelExample_mex.c'; % C ODE File
% modelfile = 'modelExample_mex.mexw64'; % MEX ODE File for Windows
% modelfile = 'modelExample_mex.mexmaci64'; % MEX ODE File for macOS
% modelfile = 'modelExample_mex.mexa64'; % MEX ODE File for Linux

% =============== Time ================ %
tspan = 0 : 0.5 : 20;

% ========= Initial Condition ========= %
y0(1) = 0; % S1
y0(2) = 0; % S2
y0(3) = 0; % S3

% ========= Parameter Values ========== %
param(1) = 0;    % S4
param(2) = 5;    % S0
param(3) = 5.5;  % J1_Vmax
param(4) = 4;    % J1_n
param(5) = 0.5;  % J1_K
param(6) = 0.1;  % J2_J2_k
param(7) = 0.1;  % J3_J3_k
param(8) = 0.01; % J0_J0_k
param(9) = 1;    % compart

% ============ ODE Solver ============= %
fast_flag = 0; % # fast_flag (0: MATLAB ODEXX)
% fast_flag = 1; % # fast_flag (1: SundialsTB CVODE)
% fast_flag = 2; % # fast_flag (2: IQMTools CVODE MEX)

% ============ Simulation ============= %
[ T, Y ] = RCGAsimulate(model, tspan, y0, param, fast_flag);

% =============== Figure ============== %
figure;
plot(T,Y,'-','LineWidth',2);
set(gca,'FontSize',10,'FontName','Arial');
legend('S1','S2','S3');
xlabel('Time');
ylabel('Concentration');
