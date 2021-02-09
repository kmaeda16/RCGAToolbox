clear all;

n_repeat = 5;


%% Start Calculation
cwd = pwd;
addpath(cwd);

for Name = {'B2','B4','B5','B6'}
    name = char(Name);
    for i = 1 : n_repeat
        dirname = sprintf('result_%s_%d',name,i);
        mkdir(dirname);
        cd(dirname);
        if strcmp('B6',name)
            copyfile('../../function/BioPreDynBenchFiles/B6/Matlab/example_optimization/dm_hkgn53_wls_5_003','.');
        end
%         Benchmark_eSS(name,i); % For normal calculation
        batch(@Benchmark_eSS,0,{name,i}); % For batch calculation (Parallel Computing Toolbox required)
        cd('..');
    end
end

rmpath(cwd);


%% Make Result files
% After finishing calculation, execute below in the directory eSS.

% dirname= 'Results';
% mkdir(dirname);
% cd(dirname);
% 
% n_repeat = 5;
% 
% 
% addpath(genpath('../../function'));
% 
% for Name = {'B2','B4','B5','B6'}
%     name = char(Name);
%     for i = 1 : n_repeat
%         fprintf('%s %d ...\n',name,i);
%         infilename = sprintf('../result_%s_%d/ess_report.mat',name,i);
%         outfilename = sprintf('eSS_%s_transition_%d.dat',name,i);
%         makeResultFiles(infilename,outfilename);
%     end
% end
% 
% rmpath(genpath('../../function'));
% 
% cd('..');
