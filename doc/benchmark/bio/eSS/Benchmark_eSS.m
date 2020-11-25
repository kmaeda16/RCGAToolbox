function Benchmark_eSS(problem_name,idum)

% idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Init
addpath(genpath('../../function'));


%% Solving
opts = [];
fprintf('\n********** %s **********\n',problem_name);
[problem, opts] = getParam(problem_name,opts);

ess_problem.f = func2str(problem.fitnessfun); %mfile containing the objective function
ess_problem.vtr = opts.vtr;
ess_problem.x_L = problem.decodingfun(zeros(1,problem.n_gene)); %lower bounds
ess_problem.x_U = problem.decodingfun( ones(1,problem.n_gene)); %upper bounds
ess_problem.c_L = -inf(1,problem.n_constraint);
ess_problem.c_U = zeros(1,problem.n_constraint);
ess_opts.maxtime = opts.maxtime;
ess_opts.maxeval = opts.maxeval;
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

opts.out_best = sprintf('eSS_%s_final_%d.dat',problem_name,idum);
writeBestAlt(elapsedTime, generation, problem, opts, x, neval);


%% Deinit
addpath(genpath('../../function'));
