% This script demonstrates how to run a real-coded genetic algorithm to
% estimate model parameters in an example kinetic model.
% 
% ------------------------ Example Kinetic Model ------------------------
% - Initial States
% S1 = 0
% S2 = 0
% S3 = 0
% 
% - Model Parameters
% S4 = 0
% S0 = 5
% J1_Vmax = 5.5
% J1_n = 4
% J1_K = 0.5
% J2_J2_k = 0.1
% J3_J3_k = 0.1
% J0_J0_k = 0.01
% compart = 1
% 
% - Reaction Kinetics
% J1 = J1_Vmax * power(S1, J1_n) / (power(J1_K, J1_n) + power(S1, J1_n))
% J2 = J2_J2_k * S2
% J3 = J3_J3_k * S3
% J0 = J0_J0_k * S0
% 
% - Differential Equations
% S1_dot = ( J0 - J1 ) / compart
% S2_dot = ( J1 - J2 ) / compart
% S3_dot = ( J2 - J3 ) / compart
% -----------------------------------------------------------------------


clearvars;

% ========= Problem Settings ========= %
modelfile = 'model_Example_odefun.m'; % Model File
decodingfun = @decoding_Example; % Decoding Function
measurement = 'measurement_Example.xls'; % Measurement File

% ========== Executing RCGA ========== %
% Results = RCGA_UNDXMGG_PE(modelfile,decodingfun,measurement); % UNDX/MGG
Results = RCGA_REXstarJGG_PE(modelfile,decodingfun,measurement); % REXstar/JGG
