clearvars;
clear randn rand;

addpath(genpath('../../../RCGA'));
addpath('../../../PE');

Param.par = 0;
Param.output_intvl = 1;
Param.vtr = 1e-6;
Param.t_limit = 1e+4;
Param.t_rextar = 6.0;
Param.selection_type = 0;
Param.out_transition = 'Transition.dat';
Param.out_population = 'Population.dat';
Param.out_solution = 'Solution.dat';

Param.n_gene = 9;
Param.n_generation = 300;
Param.n_population = 20;
Param.n_children = 20;
% Param.n_population = 1000;
% Param.n_children = 1000;
Param.output_intvl = 5;
Param.n_constraint = 0;
Param.n_parent = Param.n_gene + 1;
Param.selection_type = 0;
Param.vtr = 0;
Param.Pf = 0.45;
% Param.Pf = 0;

rng(3);

model = 'SBMLexampleLevel2.xml';
% model = SBmodel('SBMLexampleLevel2.xml');
% model = 'hill_mex';
% measurment = SBmeasurement('MeasurementExample.xls');
measurment = 'MeasurementExample.xls';
fast_flag = 1;
fitnessfun = @SSR_sbml;
Param.decodingfun = @mydecodingfun;
Param.interimreportfun = @interimreportfun_sbml;
Param.par = 0;
% optimizedmodel = REXstarJGG_sbml(model,measurment,fast_flag,Param,fitnessfun);

odefun = @hill;
icfun = @initcond;
fitnessfun = @SSR_odefun;
Param.interimreportfun = @interimreportfun_odefun;
optimizedmodel = REXstarJGG_odefun(odefun,icfun,measurment,fast_flag,Param,fitnessfun);
