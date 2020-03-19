% function doBiological_MATLAB_GA(idum)
Problem_Name = {'hiv'};
idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Init
addpath(genpath('../function'));
addpath(genpath('../../../../RCGA'));


opts = [];
fprintf('\n********** %s **********\n',char(Problem_Name));
[problem, opts] = getParam(char(Problem_Name),opts);

options = optimoptions('ga',...
    'MaxGenerations',opts.n_generation,...
    'MaxTime',opts.t_limit,...
    'FunctionTolerance',eps,...
    'MaxStallGenerations',inf,...
    'ConstraintTolerance',eps,...
    'FitnessLimit',-inf,...
    'Display','final'); % 'iter' or 'final'

%     options = optimoptions('fmincon',...
%         'MaxFunctionEvaluations',1e+5);

ObjectiveFunction  = @(x) obj_wrapper_nodecoding(problem.fitnessfun, x);
ConstraintFunction = @(x) cst_wrapper_nodecoding(problem.fitnessfun, x);

LB = problem.decodingfun(zeros(1,problem.n_gene));   % Lower bound
UB = problem.decodingfun( ones(1,problem.n_gene));   % Upper bound
tic;
[x, fval, exitflag, output] = ga(ObjectiveFunction,problem.n_gene,[],[],[],[],LB,UB,ConstraintFunction,options);
%     [x, fval, exitflag, output] = fmincon(ObjectiveFunction,rand(1,problem.n_gene),[],[],[],[],LB,UB,ConstraintFunction,options);
elapsedTime = toc;
generation = output.generations;
neval = output.funccount;
f = ObjectiveFunction(x);
g = ConstraintFunction(x);
phi = sum( max(0,g) .^2 );

fprintf('Elapsed Time = %e, f = %e, phi = %e\n',elapsedTime,f,phi);
if f < opts.vtr && phi <= eps
    fprintf('Solution found\n');
else
    fprintf('Solution NOT found\n');
end

opts.out_best = sprintf('Results/MATLAB_GA_%s_final_%d.dat',char(Problem_Name),idum);
writeBest(elapsedTime, generation, problem, opts, x, neval);


%% Deinit
rmpath(genpath('../function'));
rmpath(genpath('../../../../RCGA'));
