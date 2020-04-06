% This script was created by RCGAToolbox Mission Control PE on 03-Apr-2020

% ========= Problem Settings ========= %
modelpath = ''; % Path to Model File
addpath(modelpath);
modelfile = 'modelExample_m'; % Model File
decodingpath = ''; % Path to Decoding Function File
addpath(decodingpath);
decodingfun = @decodingExample; % Decoding Function
measurement = 'measurementExample.xls'; % Measurement File

% ========= Option Settings ========== %
opts.n_population = 50; % Population Size
opts.n_children = 25; % # Children per Generation
opts.n_parent = problem.n_gene + 1; % # Parents used for REXstar
opts.t_rexstar = 6.0; % Step-size parameter for REXstar
opts.selection_type = 0; % Parameter for JGG (0: Chosen from Children, 1: Chosen from Family)
opts.n_generation = 200; % Max # Generations
opts.t_limit = 60 * 1.000000e+00; % Max Time (sec)
opts.vtr = 0.000000e+00; % Value To Be Reached
opts.output_intvl = 1; % Output Interval Generation
opts.out_transition = 'Transition.txt'; % Transition File Name
opts.out_best = 'BestIndividual.txt'; % Best Individual File Name
opts.out_population = 'FinalPopulation.txt'; % Final Population File Name
opts.out_report = 'Report.mat'; % Report File Name
opts.n_par = 1; % # Workders
fast_flag = 0; % # fast_flag
opts.local = 0; % Local Optimizer

% ======= Setting Random Seed ======== %
rng(1); % Random Seed

% ========== Executing RCGA ========== %
clear RCGAssr;
RCGA_REXstarJGG_PE(modelfile,decodingfun,measurement,fast_flag,[],opts);

% ========== Removing Path =========== %
rmpath(modelpath);
if ~strcmp(modelpath,decodingpath)
    rmpath(decodingpath);
end
