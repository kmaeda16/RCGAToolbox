clearvars;
clear randn rand;

addpath('function');
addpath(genpath('../../../RCGA'));

Param = struct;
Param = getParam(Param,'Sphere');

Param.par = 0;
Param.output_intvl = 1;
Param.vtr = 1e-6;
Param.t_limit = 1e+4;
Param.t_rextar = 6.0;
Param.selection_type = 0;
Param.out_transition = 'Transition.dat';
Param.out_population = 'Population.dat';
Param.out_solution = 'Solution.dat';

% Param.n_generation = 1000;
Param.n_population = 500;
Param.n_children = 200;
Param.output_intvl = 20;
Param.n_constraint = 0;
Param.n_parent = Param.n_gene + 1;
Param.selection_type = 0;
Param.Pf = 0.45;
% Param.Pf = 0;

rng(3);
% [ best, Population ] = UNDXMGG(Param);
[ best, Population ] = REXstarJGG(Param);
