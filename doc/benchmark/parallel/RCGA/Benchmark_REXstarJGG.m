function Benchmark_REXstarJGG(problem_name,idum,n_par)

% idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Init
addpath(genpath('../../function'));


%% Solving
opts = [];
fprintf('\n********** %s **********\n',problem_name);
[problem, opts] = getParam(problem_name,opts);
opts.n_par = n_par;
opts.out_transition = sprintf('REXstarJGG_%s_transition_%d_%dcore.dat',problem_name,idum,n_par);
RCGA_REXstarJGG(problem,opts);


%% Deinit
rmpath(genpath('../../function'));
