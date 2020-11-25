function Benchmark_REXstarJGG(problem_name,idum)

% idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Init
addpath(genpath('../../function'));


%% Solving
opts = [];
fprintf('\n********** %s **********\n',problem_name);
[problem, opts] = getParam(problem_name,opts);
opts.out_transition = sprintf('REXstarJGG_%s_transition_%d.dat',problem_name,idum);
RCGA_REXstarJGG(problem,opts);


%% Deinit
rmpath(genpath('../../function'));
