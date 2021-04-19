% This script demonstrates how to run a real-coded genetic algorithm to
% estimate model parameters in an example kinetic model.
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

% ========= Problem Settings ========= %
modelfile = @Model_Example_odefun; % ODE file (IQM Tools format)
decodingfun = @Decoding_Example; % Decoding Function
% measurement = 'Measurement_Example.xls'; % Measurement File (EXCEL format)
measurement = 'Measurement_Example.csv'; % Measurement File (CSV)

% ========== Executing RCGA ========== %
% Results = RCGA_UNDXMGG_PE(modelfile,decodingfun,measurement); % UNDX/MGG
Results = RCGA_REXstarJGG_PE(modelfile,decodingfun,measurement); % REXstar/JGG

% ======== Convergence Curve ========= %
figure;
plot(Results.Transition.time,Results.Transition.f,'LineWidth',2);
set(gca,'FontSize',10,'FontName','Arial');
title('Convergence Curve');
xlabel('Time (sec)');
ylabel('Objective Function');
