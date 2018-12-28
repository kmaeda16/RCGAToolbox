clearvars;
clear randn rand;

addpath('function');
addpath(genpath('../../../RCGA'));

problem = struct;
opts = struct;
[ problem, opts ] = getParam(problem,opts,'g01');
% opts.t_limit = 1;

% opts.par = 0;
opts.output_intvl = 100;
% opts.vtr = 1e-6;
% opts.t_limit = 1e+4;
% opts.t_rextar = 6.0;
% opts.selection_type = 0;
% opts.out_transition = 'Transition.dat';
% opts.out_population = 'Population.dat';
% opts.out_solution = 'Solution.dat';

% opts.n_generation = 2;
% opts.n_population = 500;
% opts.n_children = 200;
% opts.output_intvl = 20;
% opts.n_constraint = 0;
% opts.n_parent = problem.n_gene + 1;
% opts.selection_type = 0;
% opts.Pf = 0.45;
% Param.Pf = 0;

rng(3);
% Results = UNDXMGG(problem,opts);
Results = REXstarJGG(problem,opts);
