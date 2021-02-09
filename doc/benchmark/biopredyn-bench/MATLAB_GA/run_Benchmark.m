clear all;

n_repeat = 5;


%% Start Calculation
cwd = pwd;
addpath(cwd);
dirname = 'Results';
mkdir(dirname);
cd(dirname);

copyfile('../../function/BioPreDynBenchFiles/B6/Matlab/example_optimization/dm_hkgn53_wls_5_003','.');

for Name = {'B2','B4','B5','B6'}
    name = char(Name);
    for i = 1 : n_repeat
%         Benchmark_MATLAB_GA(name,i); % For normal calculation
        batch(@Benchmark_MATLAB_GA,0,{name,i}); % For batch calculation (Parallel Computing Toolbox required)
    end
end

cd('..');
rmpath(cwd);


%% Save diary
% After finishing caluculation, execute below under the directory
% MATLAB_GA.

% n_repeat = 5;
% 
% jobid = 0;
% myCluster = parcluster('local');
% for Name = {'B2','B4','B5','B6'}
%     name = char(Name);
%     for j = 1 : n_repeat
%         jobid = jobid + 1;
%         jobid_str = sprintf('job%d',jobid);
%         expression = sprintf('%s = myCluster.findJob(''ID'',%d);',jobid_str,jobid);
%         eval(expression);
%         diaryname = sprintf('Results/MATLAB_GA_%s_diary_%d.txt',name,j);
%         expression = sprintf('diary(%s,''%s'')',jobid_str,diaryname);
%         eval(expression);
%     end
% end


%% Make Transition files
% After saving diaries as above, do ../convertDiary2TSV.sh on Terminal
% under the directory Results. Then, excute below in the directory
% MATLAB_GA.

% n_repeat = 5;
% 
% addpath(genpath('../function'));
% 
% for Name = {'B2','B4','B5','B6'}
%     name = char(Name);
%     for i = 1 : n_repeat
%         fprintf('%s %d ...\n',name,i);   
%         TSVfilename = sprintf('Results/MATLAB_GA_%s_TSV_%d.dat',name,i);
%         finalfilename = sprintf('Results/MATLAB_GA_%s_final_%d.dat',name,i);
%         outfilename = sprintf('Results/MATLAB_GA_%s_transition_%d.dat',name,i);
%         
%         makeResultFiles_MATLAB_GA(TSVfilename,finalfilename,outfilename);
%     end
% end
% 
% rmpath(genpath('../function'));
