function Benchmark_UNDXMGG(problem_name,idum)

% idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Init
addpath(genpath('../../function'));


%% Solving
opts = [];
fprintf('\n********** %s **********\n',problem_name);
[problem, opts] = getParam(problem_name,opts);
opts.out_transition = sprintf('UNDXMGG_%s_transition_%d.dat',problem_name,idum);
RCGA_UNDXMGG(problem,opts);


%% Deinit
addpath(genpath('../../function'));
