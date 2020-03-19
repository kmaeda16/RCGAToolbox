% function doBiological_eSS(Problem_Name,idum)
Problem_Name = {'hiv'};
idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Init
addpath(genpath('../function'));
addpath(genpath('../../../../RCGA'));


%%
opts = [];
fprintf('\n********** %s **********\n',char(Problem_Name));
[problem, opts] = getParam(char(Problem_Name),opts);

ess_problem.f = func2str(problem.fitnessfun); %mfile containing the objective function
ess_problem.vtr = opts.vtr; % f = ALLOWABLE ERROR
ess_problem.x_L = problem.decodingfun(zeros(1,problem.n_gene)); %lower bounds
ess_problem.x_U = problem.decodingfun( ones(1,problem.n_gene)); %upper bounds
ess_problem.c_L = -inf(1,problem.n_constraint);
ess_problem.c_U = zeros(1,problem.n_constraint);
ess_opts.maxtime = opts.t_limit;
ess_opts.maxeval = 1e+8;
ess_opts.tolc = 1e-30;
ess_opts.iterprint = 0;

Results = ess_kernel(ess_problem,ess_opts);
elapsedTime = Results.cpu_time;
generation = nan;
neval = Results.numeval;
x = Results.xbest;
[f, g] = problem.fitnessfun(x);
phi = sum( max(0,g) .^2 );

fprintf('Elapsed Time = %e, f = %e, phi = %e\n',elapsedTime,f,phi);
if f < opts.vtr && phi <= eps
    fprintf('Solution found\n');
else
    fprintf('Solution NOT found\n');
end

opts.out_best = sprintf('Results/eSS_%s_final_%d.dat',char(Problem_Name),idum);
writeBest(elapsedTime, generation, problem, opts, x, neval);



%% Deinit
rmpath(genpath('../function'));
rmpath(genpath('../../../../RCGA'));
