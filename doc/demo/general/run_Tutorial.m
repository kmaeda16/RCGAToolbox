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
% 
% 
% Global minimum is f = 3, g(1) = 0, g(2) = 0 at x = (-1.618, 0.6180, 0, 0,
% 0, 0, 0, 0, 0, 0) or at x = (0.6180, -1.618, 0, 0, 0, 0, 0, 0, 0, 0)
% -----------------------------------------------------------


clearvars;

% ========= Problem Settings ========= %
problem.n_gene = 10; % Number of Decision Variables
problem.n_constraint = 2; % Number of Constraints
problem.fitnessfun = @Fitness_Example; % Fitness Function
problem.decodingfun = @Decoding_Example; % Decoding Function

% ========= Option Settings ========== %
opts.n_population = 100; % Population Size
opts.n_children = 100; % Number of Children per Generation
opts.n_parent = problem.n_gene + 1; % Number of Parents Used for REXstar
opts.t_rexstar = 6.0; % Step-size Parameter for REXstar
opts.selection_type = 0; % Parameter for JGG (0: Chosen from Children, 1: Chosen from Family)
opts.Pf = 0.45; % Pf
opts.local = 0; % Local Optimizer
% opts.localopts = optimoptions(@fmincon,'ConstraintTolerance',0,'MaxFunctionEvaluations',opts.n_children,'Display','off'); % Options for Local Optimizer
opts.maxgen = 1000; % Max Number of Generations
opts.maxtime = 60; % Max Time (sec)
opts.maxeval = inf; % Max Number of fitnessfun Evaluations
opts.vtr = 0; % Value To Be Reached
opts.n_par = 1; % Number of Workers for Parallel Computation
opts.output_intvl = 10; % Output Interval Generation
opts.out_transition = 'Transition.txt'; % Transition File Name
opts.out_best = 'BestIndividual.txt'; % Best Individual File Name
opts.out_population = 'FinalPopulation.txt'; % Final Population File Name
opts.out_report = 'Report.mat'; % Report File Name
opts.interimreportfun = @RCGAdefaultinterimreportfun; % Interim Report Function
opts.finalreportfun = @RCGAdefaultfinalreportfun; % Final Report Function

% ======= Setting Random Seed ======== %
rng(0); % Random Seed

% ========== Executing RCGA ========== %
% Results = RCGA_UNDXMGG(problem,opts); % UNDX/MGG
Results = RCGA_REXstarJGG(problem,opts); % REXstar/JGG

% ======== Convergence Curve ========= %
figure;
plot(Results.Transition.time,Results.Transition.f,'LineWidth',2);
set(gca,'FontSize',10,'FontName','Arial');
title('Convergence Curve');
xlabel('Time (sec)');
ylabel('Objective Function');
