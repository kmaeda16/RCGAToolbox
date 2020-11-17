function Benchmark_RCGA(idum)
% idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Init
addpath(genpath('../function'));
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
    
    opts.Pf = 0;
    
    opts.localoptim = 0;
    opts.out_transition = sprintf('Results/UNDXMGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.localoptim,idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('Results/REXstarJGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.localoptim,idum);
    RCGA_REXstarJGG(problem,opts);

    opts.localoptim = 1;
    opts.out_transition = sprintf('Results/UNDXMGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.localoptim,idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('Results/REXstarJGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.localoptim,idum);
    RCGA_REXstarJGG(problem,opts);
    
end


%% Constrained benchmark functions
for Problem_Name = BENCHMARK2
    
    opts = [];
    fprintf('\n********** %s **********\n',char(Problem_Name));
    [problem, opts] = getParam(char(Problem_Name),opts);
    
    opts.localoptim = 0;
    
    opts.Pf = 0;
    opts.out_transition = sprintf('Results/UNDXMGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.localoptim,idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('Results/REXstarJGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.localoptim,idum);
    RCGA_REXstarJGG(problem,opts);
    
    opts.Pf = 0.45;
    opts.out_transition = sprintf('Results/UNDXMGG_%s_SR_LO%d_transition_%d.dat',char(Problem_Name),opts.localoptim,idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('Results/REXstarJGG_%s_SR_LO%d_transition_%d.dat',char(Problem_Name),opts.localoptim,idum);
    RCGA_REXstarJGG(problem,opts);
    
    opts.localoptim = 1;
    
    opts.Pf = 0;
    opts.out_transition = sprintf('Results/UNDXMGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.localoptim,idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('Results/REXstarJGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.localoptim,idum);
    RCGA_REXstarJGG(problem,opts);
    
    opts.Pf = 0.45;
    opts.out_transition = sprintf('Results/UNDXMGG_%s_SR_LO%d_transition_%d.dat',char(Problem_Name),opts.localoptim,idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('Results/REXstarJGG_%s_SR_LO%d_transition_%d.dat',char(Problem_Name),opts.localoptim,idum);
    RCGA_REXstarJGG(problem,opts);
    
end


%% Deinit
rmpath(genpath('../function'));
rmpath(genpath('../../../../RCGA'));
