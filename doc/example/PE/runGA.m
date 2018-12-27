clearvars;
clear randn rand;

addpath(genpath('../../../RCGA'));
addpath('../../../PE');

opts.par = 0;
opts.output_intvl = 1;
opts.vtr = 1e-6;
opts.t_limit = 1e+4;
opts.t_rextar = 6.0;
opts.selection_type = 0;
opts.out_transition = 'Transition.dat';
opts.out_population = 'Population.dat';
opts.out_solution = 'Solution.dat';

problem.n_gene = 9;
opts.n_generation = 300;
opts.n_population = 20;
opts.n_children = 20;
% Param.n_population = 1000;
% Param.n_children = 1000;
opts.output_intvl = 5;
problem.n_constraint = 0;
opts.n_parent = problem.n_gene + 1;
opts.selection_type = 0;
opts.vtr = 0;
opts.Pf = 0.45;
% Param.Pf = 0;

rng(3);

model = 'SBMLexampleLevel2.xml';
% model = SBmodel('SBMLexampleLevel2.xml');
% model = 'hill_mex';
% measurment = SBmeasurement('MeasurementExample.xls');
measurment = 'MeasurementExample.xls';
fast_flag = 1;
fitnessfun = @SSR_sbml;
decodingfun = @mydecodingfun;
opts.interimreportfun = @interimreportfun_sbml;
opts.par = 0;
simopts = struct;
% optimizedmodel = REXstarJGG_sbml(model,measurment,fast_flag,fitnessfun,decodingfun,simopts,opts);
% optimizedmodel = REXstarJGG_sbml(model,decodingfun,measurment);
% optimizedmodel = REXstarJGG_sbml(model,decodingfun,measurment,opts);

odefun = @hill;
icfun = @initcond;
fitnessfun = @SSR_odefun;
opts.interimreportfun = @interimreportfun_odefun;
% optimizedmodel = REXstarJGG_odefun(odefun,icfun,decodingfun,measurment,9);
optimizedmodel = REXstarJGG_odefun(odefun,icfun,decodingfun,measurment,9,opts);