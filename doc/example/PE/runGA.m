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
opts.out_best = 'Best.dat';

problem.n_gene = 9;
opts.n_generation = 50;
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

fast_flag = 0;
% fast_flag = 1;

measurment = 'MeasurementExample.xls';
% measurment = SBmeasurement('MeasurementExample.xls');
decodingfun = @mydecodingfun;

%% SBML ver
opts.method = 'ode15s';
model = 'SBMLexampleLevel2.xml';
% model = SBmodel('SBMLexampleLevel2.xml');
% optimizedmodel = REXstarJGG_sbml(model,measurment,fast_flag,fitnessfun,decodingfun,simopts,opts);
% optimizedmodel = REXstarJGG_sbml(model,decodingfun,measurment);
% optimizedmodel = REXstarJGG_sbml(model,decodingfun,measurment,[],[],[],fast_flag,[],[]);
[Results,optimizedmodel] = REXstarJGG_sbml(model,decodingfun,measurment,[],[],[],fast_flag,[],opts);

%% MATLAB ODEFUN ver
odefun = @hill;
icfun = @initcond;
% Results = REXstarJGG_odefun(odefun,icfun,decodingfun,measurment,9);
% Results = REXstarJGG_odefun(odefun,icfun,decodingfun,measurment,9,opts);
% Results = REXstarJGG_odefun(odefun,icfun,decodingfun,measurment,9,[],[],fast_flag,[],opts);

%% C ODEFUN ver
model = 'hill_c';
% Results = REXstarJGG_c(model,decodingfun,measurment,9,[],[],fast_flag,[],opts);
