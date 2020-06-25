% This script demonstrates how to run real-coded genetic algorithm to solve
% an example parameter estimation problem.


clearvars;

% ========= Problem Settings ========= %
% modelfile = 'modelExample.xml'; % SBML File
% modelfile = IQMmodel('modelExample.xml'); % Creating IQMmodel
modelfile = 'modelExample_odefun.m'; % MATLAB ODE Function File
% modelfile = 'modelExample_mex.c'; % C ODE File
% modelfile = 'modelExample_mex.mexw64'; % MEX ODE File for Windows
% modelfile = 'modelExample_mex.mexmaci64'; % MEX ODE File for macOS
% modelfile = 'modelExample_mex.mexa64'; % MEX ODE File for Linux
decodingfun = @decodingExample; % Decoding Function
measurement = 'measurementExample.xls'; % Measurement File

% ========= Option Settings ========== %
opts.n_population = 50; % Population Size
opts.n_children = 25; % # Children per Generation
% opts.n_parent = n_gene + 1; % # Parents used for REXstar
opts.t_rexstar = 6.0; % Step-size parameter for REXstar
opts.selection_type = 0; % Parameter for JGG (0: Chosen from Children, 1: Chosen from Family)
opts.n_generation = 200; % Max # Generations
opts.maxtime = 60 * 1.000000e+00; % Max Time (sec)
opts.maxeval = inf; %% Max # fitnessfun Evaluations
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
