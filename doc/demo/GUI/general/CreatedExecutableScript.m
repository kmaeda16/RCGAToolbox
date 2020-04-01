% ========= Problem Settings ========= %
problem.n_gene = 10; % # Variables
problem.n_constraint = 2; % # Constraints
fitnesspath = ''; % Path to Fitness Function File
addpath(fitnesspath);
problem.fitnessfun = @fitnessExample; % Fitness Function
decodingpath = ''; % Path to Decoding Function File
addpath(decodingpath);
problem.decodingfun = @decodingExample; % Decoding Function

% ========= Option Settings ========== %
opts.n_population = 200; % Population Size
opts.n_children = 100; % # Children per Generation
opts.n_parent = problem.n_gene + 1; % # Parents used for REXstar
opts.t_rexstar = 6.0; % Step-size parameter for REXstar
opts.selection_type = 0; % Parameter for JGG (0: Chosen from Children, 1: Chosen from Family)
opts.Pf = 4.500000e-01; % Pf
opts.n_generation = 1000; % Max # Generations
opts.t_limit = 60 * 1.000000e+00; % Max Time (sec)
opts.vtr = 0.000000e+00; % Value To Be Reached
opts.output_intvl = 10; % Output Interval Generation
opts.out_transition = 'Transition.txt'; % Transition File Name
opts.out_best = 'BestIndividual.txt'; % Best Individual File Name
opts.out_population = 'FinalPopulation.txt'; % Final Population File Name
opts.out_report = 'Report.mat'; % Report File Name
opts.n_par = 1; % # Workders
opts.local = 0; % Local Optimizer

% ======= Setting Random Seed ======== %
rng(1); % Random Seed

% ========== Executing RCGA ========== %
RCGA_REXstarJGG(problem,opts);

% ========== Removing Path =========== %
rmpath(fitnesspath);
if ~strcmp(fitnesspath,decodingpath)
    rmpath(decodingpath);
end
