function doBiological_MATLAB_GA(problem_name,idum)
% problem_name = 'hiv';
% idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Init
% addpath(genpath('../function'));
% addpath(genpath('../../../../RCGA'));
addpath(genpath('../../function'));
addpath(genpath('../../../../../RCGA'));

%%
opts = [];

fprintf('\n********** %s **********\n',problem_name);
[problem, opts] = getParam(problem_name,opts);
opts.out_best = sprintf('MATLAB_GA_%s_final_%d.dat',problem_name,idum);

options = optimoptions('ga',...
    'MaxGenerations',opts.n_generation,...
    'MaxTime',opts.t_limit,...
    'FunctionTolerance',eps,...
    'MaxStallGenerations',inf,...
    'ConstraintTolerance',eps,...
    'FitnessLimit',-inf,...
    'Display','iter'); % 'iter' or 'final'

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

opts.out_best = sprintf('MATLAB_GA_%s_final_%d.dat',problem_name,idum);
writeBest(elapsedTime, generation, problem, opts, x, neval);


%% Deinit
% rmpath(genpath('../function'));
% rmpath(genpath('../../../../RCGA'));
addpath(genpath('../../function'));
addpath(genpath('../../../../../RCGA'));

