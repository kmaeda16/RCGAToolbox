function doBenchmark_RCGA(idum)


%% Init
addpath(genpath('function'));
addpath(genpath('../../../RCGA'));

BENCHMARK1 = {'Sphere','ScaledSphere','Ellipsoid','Cigar','k_tablet', ...
    'MMbenchmark','Rosenbrock_star','Rosenbrock_chain','Ackley','Bohachevsky', ...
    'Rastrigin','Schaffer','Schwefel'};
BENCHMARK2 = {'g01','g02','g03','g04','g05', ...
    'g06','g07', 'g08','g09','g10', ...
    'g11','g12','g13'};

rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Unconstrained benchmark functions
opts.Pf = 0;
for Problem_Name = BENCHMARK1
    
    fprintf('\n********** %s **********\n',char(Problem_Name));
    [problem, opts] = getParam(char(Problem_Name),opts);
%     opts.vtr             = 5e+2;
    
    opts.out_transition = sprintf('UNDXMGG_%s_NR_transition_%d.dat',char(Problem_Name),idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('REXstarJGG_%s_NR_transition_%d.dat',char(Problem_Name),idum);
    RCGA_REXstarJGG(problem,opts);
    
end


%% Constrained benchmark functions
for Problem_Name = BENCHMARK2
    opts = [];
    fprintf('\n********** %s **********\n',char(Problem_Name));
    [problem, opts] = getParam(char(Problem_Name),opts);
    
    opts.Pf = 0;
    opts.out_transition = sprintf('UNDXMGG_%s_NR_transition_%d.dat',char(Problem_Name),idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('REXstarJGG_%s_NR_transition_%d.dat',char(Problem_Name),idum);
    RCGA_REXstarJGG(problem,opts);
    
    opts.Pf = 0.45;
    opts.out_transition = sprintf('UNDXMGG_%s_SR_transition_%d.dat',char(Problem_Name),idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('REXstarJGG_%s_SR_transition_%d.dat',char(Problem_Name),idum);
    RCGA_REXstarJGG(problem,opts);
    
end


%% Deinit
rmpath(genpath('function'));
rmpath(genpath('../../../RCGA'));
