% This script demonstrates how to run a real-coded genetic algorithm to
% solve an example constrained optimization problem.
% 
% --------------------- Example Problem ---------------------
% Minimize:
%   f = x(1)^2 + x(2)^2 + ... + x(10)^2
% 	
% Subject to:
%   g(1) = x(1) * x(2) + 1 <= 0
%   g(2) = x(1) + x(2) + 1 <= 0
% 	-5.12 <= x(i) <= 5.12 for all i
% -----------------------------------------------------------


clearvars;

% ========= Problem Settings ========= %
problem.n_gene = 10; % # Decision Variables
problem.n_constraint = 2; % # Constraints
problem.fitnessfun = @fitnessExample; % Fitness Function
problem.decodingfun = @decodingExample; % Decoding Function

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
