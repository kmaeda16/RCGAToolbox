% function doBiological_RCGA(Problem_Name,idum)
Problem_Name = {'hiv'};
idum = 1;
rng(idum); % For Reproducibility
fprintf('idum = %d\n',idum);

%% Init
addpath(genpath('../function'));
addpath(genpath('../../../../RCGA'));


%%
opts = [];
fprintf('\n********** %s **********\n',char(Problem_Name));
[problem, opts] = getParam(char(Problem_Name),opts);
opts.out_transition = sprintf('Results/REXstarJGG_%s_transition_%d.dat',char(Problem_Name),idum);
% RCGA_REXstarJGG(problem,opts);


%% Deinit
rmpath(genpath('../function'));
rmpath(genpath('../../../../RCGA'));
