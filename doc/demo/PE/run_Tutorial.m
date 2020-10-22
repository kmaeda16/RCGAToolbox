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
% modelfile = 'model_Example.xml'; % SBML File
% modelfile = IQMmodel('model_Example.xml'); % Creating IQMmodel
modelfile = 'model_Example_odefun.m'; % MATLAB ODE Function File
% modelfile = 'model_Example_mex.c'; % C ODE File
% modelfile = 'model_Example_mex.mexw64'; % MEX ODE File for Windows
% modelfile = 'model_Example_mex.mexmaci64'; % MEX ODE File for macOS
% modelfile = 'model_Example_mex.mexa64'; % MEX ODE File for Linux
decodingfun = @decoding_Example; % Decoding Function
measurement = 'measurement_Example.xls'; % Measurement File

% ========= Option Settings ========== %
opts.n_population = 50; % Population Size
opts.n_children = 25; % # Children per Generation
% opts.n_parent = n_gene + 1; % # Parents used for REXstar
opts.t_rexstar = 6.0; % Step-size parameter for REXstar
opts.selection_type = 0; % Parameter for JGG (0: Chosen from Children, 1: Chosen from Family)
opts.n_generation = 200; % Max # Generations
opts.maxtime = 60 * 1.000000e+00; % Max Time (sec)
opts.maxeval = inf; % Max # fitnessfun Evaluations
opts.vtr = 0.000000e+00; % Value To Be Reached
opts.output_intvl = 1; % Output Interval Generation
opts.out_transition = 'Transition.txt'; % Transition File Name
opts.out_best = 'BestIndividual.txt'; % Best Individual File Name
opts.out_population = 'FinalPopulation.txt'; % Final Population File Name
opts.out_report = 'Report.mat'; % Report File Name
opts.n_par = 1; % # Workders
fast_flag = 0; % # fast_flag (0: MATLAB ODEXX)
% fast_flag = 1; % # fast_flag (1: SundialsTB CVODE)
% fast_flag = 2; % # fast_flag (2: IQMTools CVODE MEX)
opts.local = 0; % Local Optimizer

% ======= Setting Random Seed ======== %
rng(1); % Random Seed

% ========== Executing RCGA ========== %
clear RCGAssr;
% Results = RCGA_UNDXMGG_PE(modelfile,decodingfun,measurement,fast_flag,[],opts); % UNDX/MGG
Results = RCGA_REXstarJGG_PE(modelfile,decodingfun,measurement,fast_flag,[],opts); % REXstar/JGG
