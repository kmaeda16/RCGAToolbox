function doBiological_RCGA_REXstarJGG(problem_name,idum)
% problem_name = 'hiv';
% idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);


%% Init
% addpath(genpath('../function'));
% addpath(genpath('../../../../RCGA'));
addpath(genpath('../../function'));
addpath(genpath('../../../../../RCGA'));


%%
opts = [];
fprintf('\n********** %s **********\n',problem_name);
[problem, opts] = getParam(problem_name,opts);
opts.out_transition = sprintf('REXstarJGG_%s_transition_%d.dat',problem_name,idum);
RCGA_REXstarJGG(problem,opts);


%% Deinit
% rmpath(genpath('../function'));
% rmpath(genpath('../../../../RCGA'));
addpath(genpath('../../function'));
addpath(genpath('../../../../../RCGA'));
