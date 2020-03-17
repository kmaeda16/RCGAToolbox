% function doBenchmark_RCGA(idum)
idum = 1;


%% Init
addpath(genpath('../function'));
addpath(genpath('../../../../RCGA'));

BENCHMARK1 = {    'MMbenchmark'};
BENCHMARK2 = {'g13'};

rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Unconstrained benchmark functions
opts.Pf = 0;
for Problem_Name = BENCHMARK1
    
    fprintf('\n********** %s **********\n',char(Problem_Name));
    [problem, opts] = getParam(char(Problem_Name),opts);
%     opts.vtr             = 5e+2;
     opts.n_localoptimind = 1;
%     opts.out_transition = sprintf('UNDXMGG_%s_NR_transition_%d.dat',char(Problem_Name),idum);
%     RCGA_UNDXMGG(problem,opts);
%     opts.out_transition = sprintf('REXstarJGG_%s_NR_transition_%d.dat',char(Problem_Name),idum);
%     RCGA_REXstarJGG(problem,opts);
    
end


%% Constrained benchmark functions
for Problem_Name = BENCHMARK2
    opts = [];
    fprintf('\n********** %s **********\n',char(Problem_Name));
    [problem, opts] = getParam(char(Problem_Name),opts);
    
    opts.Pf = 0;
%     opts.out_transition = sprintf('UNDXMGG_%s_NR_transition_%d.dat',char(Problem_Name),idum);
%     RCGA_UNDXMGG(problem,opts);
%     opts.out_transition = sprintf('REXstarJGG_%s_NR_transition_%d.dat',char(Problem_Name),idum);
%     RCGA_REXstarJGG(problem,opts);
    
    opts.Pf = 0.45;
    opts.n_localoptimind = 1;
%     opts.out_transition = sprintf('UNDXMGG_%s_SR_transition_%d.dat',char(Problem_Name),idum);
%     RCGA_UNDXMGG(problem,opts);
    opts.out_transition = sprintf('REXstarJGG_%s_SR_transition_%d.dat',char(Problem_Name),idum);
    RCGA_REXstarJGG(problem,opts);
    
end


%% Deinit
rmpath(genpath('../function'));
rmpath(genpath('../../../../RCGA'));
