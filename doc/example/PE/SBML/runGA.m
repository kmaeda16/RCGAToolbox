clearvars;
clear randn rand;

addpath(genpath('../../../../RCGA'));
addpath('../../../../PE');
addpath('../Common');

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

% fast_flag = 0; % ODEXX
% fast_flag = 1; % SundialsTB
fast_flag = 2; % MEX

measurement = 'MeasurementExample.xls';
% measurment = SBmeasurement('MeasurementExample.xls');
decodingfun = @mydecodingfun;

%% SBML ver
% simopts.method = 'ode23s';
simopts.abstol = 1e-6;
simopts.reltol = 1e-6;
model = 'SBMLexampleLevel2.xml';
model = 'SBMLexampleLevel2';
% model = IQMmodel('SBMLexampleLevel2.xml');
% optimizedmodel = REXstarJGG_PE(model,measurment,fast_flag,fitnessfun,decodingfun,simopts,opts);
% optimizedmodel = REXstarJGG_PE(model,decodingfun,measurement);
% optimizedmodel = REXstarJGG_PE(model,decodingfun,measurement,[],[],[],fast_flag,[],[]);
Results = REXstarJGG_PE(model,decodingfun,measurement,[],[],fast_flag,simopts,opts);

%% MATLAB ODEFUN ver
% odefun = @hill_odefun;
% Results = REXstarJGG_PE(odefun,icfun,decodingfun,measurment,9);
% Results = REXstarJGG_PE(odefun,icfun,decodingfun,measurment,9,opts);
% Results = REXstarJGG_PE(odefun,decodingfun,measurement,[],[],fast_flag,[],opts);

%% C ODEFUN ver
% model = 'hill_c';
% Results = REXstarJGG_PE(model,decodingfun,measurement,[],[],fast_flag,[],opts);
