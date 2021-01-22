function Benchmark_eSS(idum)

% idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Init
addpath(genpath('../function'));

BENCHMARK1 = {'Sphere','ScaledSphere','Ellipsoid','Cigar','k_tablet', ...
    'MMbenchmark','Rosenbrock_star','Rosenbrock_chain','Ackley','Bohachevsky', ...
    'Rastrigin','Schaffer','Schwefel','ScaledShiftedRotatedRastrigin'};
BENCHMARK2 = {'g01','g02','g03','g04','g05', ...
    'g06','g07', 'g08','g09','g10', ...
    'g11','g12','g13'};


%% Unconstrained problems
for Problem_Name = BENCHMARK1
    
    opts = [];
    fprintf('\n********** %s **********\n',char(Problem_Name));
    [problem, opts] = getParam(char(Problem_Name),opts);
    
    ess_problem = [];
    ess_opts = [];
    ess_problem.f = char(Problem_Name); %mfile containing the objective function
    ess_problem.vtr = opts.vtr;
    ess_problem.x_L = problem.decodingfun(zeros(1,problem.n_gene)); %lower bounds
    ess_problem.x_U = problem.decodingfun( ones(1,problem.n_gene)); %upper bounds
    ess_opts.maxtime = opts.maxtime;
    ess_opts.maxeval = opts.maxeval;
    ess_opts.tolc = 1e-30;
    ess_opts.iterprint = 0;
    ess_opts.inter_save = 0;
    ess_opts.local.solver = 'fmincon';
    ess_opts.local.finish = 'fmincon';
    
    Results = ess_kernel(ess_problem,ess_opts);
    
    elapsedTime = Results.cpu_time;
    generation = nan;
    neval = Results.numeval;
    x = Results.xbest;
    f = Results.fbest;
    
    fprintf('Elapsed Time = %e, f = %e\n',elapsedTime,f);
    if f < opts.vtr
        fprintf('Solution found\n');
    else
        fprintf('Solution NOT found\n');
    end
    
    opts.out_best = sprintf('Results/eSS_%s_final_%d.dat',char(Problem_Name),idum);
    writeBestAlt(elapsedTime, generation, problem, opts, x, neval);
    
end

%% Constrained problems
for Problem_Name = BENCHMARK2
    
    opts = [];
    fprintf('\n********** %s **********\n',char(Problem_Name));
    [problem, opts] = getParam(char(Problem_Name),opts);
    
    ess_problem = [];
    ess_opts = [];
    ess_problem.f = char(Problem_Name); %mfile containing the objective function
    ess_problem.vtr = opts.vtr;
    ess_problem.x_L = problem.decodingfun(zeros(1,problem.n_gene)); %lower bounds
    ess_problem.x_U = problem.decodingfun( ones(1,problem.n_gene)); %upper bounds
    ess_problem.c_L = -inf(1,problem.n_constraint);
    ess_problem.c_U = zeros(1,problem.n_constraint);
    ess_opts.maxtime = opts.maxtime;
    ess_opts.maxeval = opts.maxeval;
    ess_opts.tolc = 1e-30;
    ess_opts.iterprint = 0;
    ess_opts.inter_save = 0;
    ess_opts.local.solver = 'fmincon';
    ess_opts.local.finish = 'fmincon';
    
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
    writeBestAlt(elapsedTime, generation, problem, opts, x, neval);
    
end


%% Deinit
rmpath(genpath('../function'));
