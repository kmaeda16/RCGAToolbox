function Benchmark_MATLAB_GA(idum)
% idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Init
addpath(genpath('../function'))
addpath(genpath('../../../../RCGA'));

BENCHMARK1 = {'Sphere','ScaledSphere','Ellipsoid','Cigar','k_tablet', ...
    'MMbenchmark','Rosenbrock_star','Rosenbrock_chain','Ackley','Bohachevsky', ...
    'Rastrigin','Schaffer','Schwefel'};
BENCHMARK2 = {'g01','g02','g03','g04','g05', ...
    'g06','g07', 'g08','g09','g10', ...
    'g11','g12','g13'};


%% Unconstrained benchmark functions
for Problem_Name = BENCHMARK1
    
    opts = [];
    fprintf('\n********** %s **********\n',char(Problem_Name));
    [problem, opts] = getParam(char(Problem_Name),opts);
    
    options = optimoptions('ga',...
        'MaxGenerations',opts.n_generation,...
        'MaxTime',opts.t_limit,...
        'FunctionTolerance',eps,...
        'MaxStallGenerations',inf,...
        'ConstraintTolerance',eps,...
        'FitnessLimit',opts.vtr,...
        'Display','final'); % 'iter' or 'final'

%     options = optimoptions('fmincon',...
%         'MaxFunctionEvaluations',1e+5);
        
    ObjectiveFunction = @(x) obj_wrapper_nodecoding(problem.fitnessfun, x);
    
    LB = problem.decodingfun(zeros(1,problem.n_gene));   % Lower bound
    UB = problem.decodingfun( ones(1,problem.n_gene));   % Upper bound
    tic;
    [x, fval, exitflag, output] = ga(ObjectiveFunction,problem.n_gene,[],[],[],[],LB,UB,[],options);
%     [x, fval, exitflag, output] = fmincon(ObjectiveFunction,rand(1,problem.n_gene),[],[],[],[],LB,UB,[],options); output.generations = nan;
    elapsedTime = toc;
    generation = output.generations;
    neval = output.funccount;
    f = ObjectiveFunction(x);
    
    fprintf('Elapsed Time = %e, f = %e\n',elapsedTime,f);
    if f < opts.vtr
        fprintf('Solution found\n');
    else
        fprintf('Solution NOT found\n');
    end
    
    opts.out_best = sprintf('Results/MATLAB_GA_%s_final_%d.dat',char(Problem_Name),idum);
    writeBest(elapsedTime, output.generations, problem, opts, x, neval);
    
end


%% Constrained benchmark functions
for Problem_Name = BENCHMARK2
    
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

end


%% Deinit
rmpath(genpath('../function'));
rmpath(genpath('../../../../RCGA'));
