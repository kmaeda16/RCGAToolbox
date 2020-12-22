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


clearvars;

% ========= Problem Settings ========= %
% modelfile = 'Model_Example_SBML.xml'; % SBML file (IQM Tools required)
% modelfile = IQMmodel('Model_Example_SBML.xml'); % Creating an IQMmodel object (IQM Tools required)
modelfile = @Model_Example_odefun; % ODE file (IQM Tools format)
% modelfile = 'Model_Example_odefun.m'; % ODE file (IQM Tools format)
% modelfile = 'Model_Example_SBML_mex.c'; % C source code (IQM Tools required)
% modelfile = 'Model_Example_SBML_mex.mexw64'; % MEX model for Windows
% modelfile = 'Model_Example_SBML_mex.mexmaci64'; % MEX model file for macOS
% modelfile = 'Model_Example_SBML_mex.mexa64'; % MEX model file for Linux
decodingfun = @Decoding_Example; % Decoding Function
% measurement = 'Measurement_Example.xls'; % Measurement File (EXCEL format)
measurement = 'Measurement_Example.csv'; % Measurement File (CSV)

% ========= Option Settings ========== %
opts.n_population = 50; % Population Size
opts.n_children = 25; % # Children per Generation
% opts.n_parent = n_gene + 1; % # Parents Used for REXstar
opts.t_rexstar = 6.0; % Step-size parameter for REXstar
opts.selection_type = 0; % Parameter for JGG (0: Chosen from Children, 1: Chosen from Family)
opts.maxgen = 200; % Max # Generations
opts.maxtime = 60; % Max Time (sec)
opts.maxeval = inf; % Max # fitnessfun Evaluations
opts.vtr = 0.000000e+00; % Value To Be Reached
opts.output_intvl = 1; % Output Interval Generation
opts.out_transition = 'Transition.txt'; % Transition File Name
opts.out_best = 'BestIndividual.txt'; % Best Individual File Name
opts.out_population = 'FinalPopulation.txt'; % Final Population File Name
opts.out_report = 'Report.mat'; % Report File Name
opts.n_par = 1; % # Workers for Parallel Computation
fast_flag = 0; % # fast_flag (0: MATLAB ODEXX)
% fast_flag = 1; % # fast_flag (1: SundialsTB CVODE) (SundialsTB required)
% fast_flag = 2; % # fast_flag (2: IQM Tools CVODE MEX) (IQM Tools required)
opts.local = 0; % Local Optimizer

% ======= Setting Random Seed ======== %
rng(0); % Random Seed

% ========== Executing RCGA ========== %
clear RCGAssr;
% Results = RCGA_UNDXMGG_PE(modelfile,decodingfun,measurement,fast_flag,[],opts); % UNDX/MGG
Results = RCGA_REXstarJGG_PE(modelfile,decodingfun,measurement,fast_flag,[],opts); % REXstar/JGG
