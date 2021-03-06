function Benchmark_RCGA(idum)

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
    
    opts.Pf = 0;
    
    opts.local = 0;
    opts.out_transition = sprintf('Results/UNDXMGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.local,idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('Results/REXstarJGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.local,idum);
    RCGA_REXstarJGG(problem,opts);

    opts.local = 1;
    opts.out_transition = sprintf('Results/UNDXMGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.local,idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('Results/REXstarJGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.local,idum);
    RCGA_REXstarJGG(problem,opts);
    
end


%% Constrained problems
for Problem_Name = BENCHMARK2
    
    opts = [];
    fprintf('\n********** %s **********\n',char(Problem_Name));
    [problem, opts] = getParam(char(Problem_Name),opts);
    
    opts.local = 0;
    
    opts.Pf = 0;
    opts.out_transition = sprintf('Results/UNDXMGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.local,idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('Results/REXstarJGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.local,idum);
    RCGA_REXstarJGG(problem,opts);
    
    opts.Pf = 0.45;
    opts.out_transition = sprintf('Results/UNDXMGG_%s_SR_LO%d_transition_%d.dat',char(Problem_Name),opts.local,idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('Results/REXstarJGG_%s_SR_LO%d_transition_%d.dat',char(Problem_Name),opts.local,idum);
    RCGA_REXstarJGG(problem,opts);
    
    opts.local = 1;
    
    opts.Pf = 0;
    opts.out_transition = sprintf('Results/UNDXMGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.local,idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('Results/REXstarJGG_%s_NR_LO%d_transition_%d.dat',char(Problem_Name),opts.local,idum);
    RCGA_REXstarJGG(problem,opts);
    
    opts.Pf = 0.45;
    opts.out_transition = sprintf('Results/UNDXMGG_%s_SR_LO%d_transition_%d.dat',char(Problem_Name),opts.local,idum);
    RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('Results/REXstarJGG_%s_SR_LO%d_transition_%d.dat',char(Problem_Name),opts.local,idum);
    RCGA_REXstarJGG(problem,opts);
    
end


%% Deinit
rmpath(genpath('../function'));
