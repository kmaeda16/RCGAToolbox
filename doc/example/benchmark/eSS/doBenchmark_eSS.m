% function doBenchmark_eSS(idum)
idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Init
addpath(genpath('../function'));
addpath(genpath('../../../../RCGA'));

BENCHMARK1 = {'Sphere','ScaledSphere','Ellipsoid','Cigar','k_tablet', ...
    'MMbenchmark','Rosenbrock_star','Rosenbrock_chain','Ackley','Bohachevsky', ...
    'Rastrigin','Schaffer','Schwefel'};
% BENCHMARK1 = {'Rosenbrock_chain'};
BENCHMARK2 = {'g01','g02','g03','g04','g05', ...
    'g06','g07', 'g08','g09','g10', ...
    'g11','g12','g13'};
% BENCHMARK2 = {'g02'};


opts = [];


%% Unconstrained benchmark functions
% for Problem_Name = BENCHMARK1
%     
%     fprintf('\n********** %s **********\n',char(Problem_Name));
% 
%     opts = [];
%     [problem, opts] = getParam(char(Problem_Name),opts);
%     ess_problem.f = char(Problem_Name); %mfile containing the objective function
%     ess_problem.vtr = opts.vtr; % f = ALLOWABLE ERROR
%     
%     ess_problem.x_L = problem.decodingfun(zeros(1,problem.n_gene)); %lower bounds
%     ess_problem.x_U = problem.decodingfun( ones(1,problem.n_gene)); %upper bounds
%     ess_problem.c_L = -inf(1,problem.n_constraint);
%     ess_problem.c_U = zeros(1,problem.n_constraint);
%     
%     ess_opts.maxtime = opts.t_limit;
%     ess_opts.inter_save = 1;
%     ess_opts.maxeval = 1e+8;
%     ess_opts.tolc = 1e-30;
%     
%     tic;
%     Results = ess_kernel(ess_problem,ess_opts);
%     elapsedTime = toc;
%     
%     fprintf('Elapsed Time = %e, f = %e\n',elapsedTime,Results.fbest);
%     if Results.fbest < opts.vtr
%         fprintf('Solution found\n');
%     else
%         fprintf('Solution NOT found\n');
%     end
%     
%     opts.out_best = sprintf('Results/eSS_%s_final_%d.dat',char(Problem_Name),idum);
%     writeBest(elapsedTime, nan, problem, opts, Results.xbest);
%     
% end

%% Constrained benchmark functions
for Problem_Name = BENCHMARK2
    
    fprintf('\n********** %s **********\n',char(Problem_Name));

    opts = [];
    [problem, opts] = getParam(char(Problem_Name),opts);
    ess_problem.f = char(Problem_Name); %mfile containing the objective function
    ess_problem.vtr = opts.vtr; % f = ALLOWABLE ERROR
    
    ess_problem.x_L = problem.decodingfun(zeros(1,problem.n_gene)); %lower bounds
    ess_problem.x_U = problem.decodingfun( ones(1,problem.n_gene)); %upper bounds
    ess_problem.c_L = -inf(1,problem.n_constraint);
    ess_problem.c_U = zeros(1,problem.n_constraint);
    
    ess_opts.maxtime = opts.t_limit;
    ess_opts.inter_save = 1;
    ess_opts.maxeval = 1e+8;
    ess_opts.tolc = 1e-30;
    ess_opts.iterprint = 0;
    
    tic;
    Results = ess_kernel(ess_problem,ess_opts);
    elapsedTime = toc;
    
    [f, g] = problem.fitnessfun(Results.xbest);
    phi = sum( max(0,g) .^2 );
    fprintf('Elapsed Time = %e, f = %e, phi = %e\n',elapsedTime,f,phi);
    if f < opts.vtr && phi <= eps
        fprintf('Solution found\n');
    else
        fprintf('Solution NOT found\n');
    end
    
    opts.out_best = sprintf('Results/eSS_%s_final_%d.dat',char(Problem_Name),idum);
    writeBest(elapsedTime, nan, problem, opts, Results.xbest);
    
end


%% Deinit
rmpath(genpath('../function'));
rmpath(genpath('../../../../RCGA'));
