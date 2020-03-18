% function main_hiv(idum)
addpath(genpath('../function'))
addpath(genpath('../../../../RCGA'));

% clear mex;
% clear all;
% close all;
	
%%
% addpath('expdata');
% 
% global ExpData;
% 
% for i = 1:16
%     filename = strcat('pseudoexpdata',num2str(i,'%d'),'.dat');
%     ExpData(i,:,:) = dlmread(filename);
% end
idum = 1;
rng(idum);

problem = [];



%========================= PROBLEM SPECIFICATIONS ===========================
problem.f = 'Sphere'; %mfile containing the objective function
decodingfun = @Sphere_decode;
n_gene = 100;

problem.f = 'ScaledSphere'; %mfile containing the objective function
decodingfun = @ScaledSphere_decode;
n_gene = 100;

problem.f = 'k_tablet'; %mfile containing the objective function
decodingfun = @k_tablet_decode;
n_gene = 100;

% problem.f = 'MMbenchmark'; %mfile containing the objective function
% decodingfun = @MMbenchmark_decode;
% n_gene = 100;

% problem.f = 'Rosenbrock_star'; %mfile containing the objective function
% decodingfun = @Rosenbrock_star_decode;
% n_gene = 50;

problem.f = 'Rosenbrock_chain'; %mfile containing the objective function
decodingfun = @Rosenbrock_chain_decode;
n_gene = 100;
n_costraint = 0;
vtr = 1e-6;

% problem.f = 'Ackley'; %mfile containing the objective function
% decodingfun = @Ackley_decode;
% n_gene = 50;

% problem.f = 'Bohachevsky'; %mfile containing the objective function
% decodingfun = @Bohachevsky_decode;
% n_gene = 50;

% 
% problem.f = 'Rastrigin'; %mfile containing the objective function
% decodingfun = @Rastrigin_decode;
% n_gene = 20;

problem.f = 'Schaffer'; %mfile containing the objective function
decodingfun = @Schaffer_decode;
n_gene = 20;
n_costraint = 0;
vtr = 1e-6;

% problem.f = 'Schwefel'; %mfile containing the objective function
% decodingfun = @Schwefel_decode;
% n_gene = 10;

% problem.f = 'g01'; %mfile containing the objective function
% decodingfun = @g01_decode;
% n_gene = 13;
% n_costraint = 9;

% problem.f = 'g02'; %mfile containing the objective function
% decodingfun = @g02_decode;
% n_gene = 20;
% n_costraint = 2;
% vtr = -0.803619;


% problem.f = 'g03'; %mfile containing the objective function
% decodingfun = @g03_decode;
% n_gene = 10;
% n_costraint = 1;

% problem.f = 'g05'; %mfile containing the objective function
% decodingfun = @g05_decode;
% n_gene = 4;
% n_costraint = 5;

% problem.f = 'g10'; %mfile containing the objective function
% decodingfun = @g10_decode;
% n_gene = 8;
% n_costraint = 6;

% problem.f = 'g13'; %mfile containing the objective function
% decodingfun = @g13_decode;
% n_gene = 5;
% n_costraint = 3;

                                                                


problem.vtr = vtr; % f = ALLOWABLE ERROR



problem.x_L = decodingfun(zeros(1,n_gene)); %lower bounds
problem.x_U = decodingfun( ones(1,n_gene)); %upper bounds

for i = 1:n_costraint
    problem.c_L(1,i) = -inf;
    problem.c_U(1,i) = 0;
end

opts = [];
opts.maxtime = 7 * 24 * 60 * 60;
opts.inter_save = 1;
opts.maxeval = 1e+8;
% opts.ndiverse = 40;
% opts.maxtime=7;
% opts.local.solver='ipopt';
% opts.local.solver='misqpg';
opts.local.solver = 'fmincon';
% opts.local.solver = 'dhc';
% opts.local.finish = 'fmincon';
% opts.local.iterprint = 1;
opts.tolc = 1e-30;
% opts.local.MaxFunctionEvaluations = 5000;
%========================= END OF PROBLEM SPECIFICATIONS =====================

Results = ess_kernel(problem,opts);

% [f, g] = wrapper_hiv_con_mex(Results.xbest);
% disp(f);
% disp(g);

% filename = sprintf('ess_report_%d.mat',idum);
% save(filename);
