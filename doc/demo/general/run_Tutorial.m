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
% 
% - Reaction Kinetics
% J1 = J1_Vmax * power(S1, J1_n) / (power(J1_K, J1_n) + power(S1, J1_n))
% J2 = J2_J2_k * S2
% J3 = J3_J3_k * S3
% J0 = J0_J0_k * S0
% 
% - Differential Equations
% S1_dot = J0 - J1
% S2_dot = J1 - J2
% S3_dot = J2 - J3
% -----------------------------------------------------------------------


clearvars;

% ========= Problem Settings ========= %
problem.n_gene = 10; % # Decision Variables
problem.n_constraint = 2; % # Constraints
problem.fitnessfun = @fitness_Example; % Fitness Function
problem.decodingfun = @decoding_Example; % Decoding Function

% ========= Option Settings ========== %
opts.n_population = 200; % Population Size
opts.n_children = 100; % # Children per Generation
opts.n_parent = problem.n_gene + 1; % # Parents used for REXstar
opts.t_rexstar = 6.0; % Step-size parameter for REXstar
opts.selection_type = 0; % Parameter for JGG (0: Chosen from Children, 1: Chosen from Family)
opts.Pf = 4.500000e-01; % Pf
opts.n_generation = 1000; % Max # Generations
opts.maxtime = 60 * 1.000000e+00; % Max Time (sec)
opts.maxeval = inf; % Max # fitnessfun Evaluations
opts.vtr = 0.000000e+00; % Value To Be Reached
opts.output_intvl = 10; % Output Interval Generation
opts.out_transition = 'Transition.txt'; % Transition File Name
opts.out_best = 'BestIndividual.txt'; % Best Individual File Name
opts.out_population = 'FinalPopulation.txt'; % Final Population File Name
opts.out_report = 'Report.mat'; % Report File Name
opts.n_par = 1; % # Workers
opts.local = 0; % Local Optimizer

% ======= Setting Random Seed ======== %
rng(1); % Random Seed

% ========== Executing RCGA ========== %
% Results = RCGA_UNDXMGG(problem,opts); % UNDX/MGG
Results = RCGA_REXstarJGG(problem,opts); % REXstar/JGG
