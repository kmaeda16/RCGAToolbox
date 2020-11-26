
n_repeat = 5;
dirname = 'Results';

mkdir(dirname);
cd(dirname);
addpath('..');


%% Start Calculation
for Name = {'hiv','threestep'}
    name = char(Name);
    for i = 1 : n_repeat
%         Benchmark_MATLAB_GA(name,i); % For normal calculation
        batch(@Benchmark_MATLAB_GA,0,{name,i}); % For batch calculation
    end
end

rmpath('..');
cd('..');


%% Save diary
% After finishing caluculation, execute below under the directory
% MATLAB_GA.
% 
% jobid = 0;
% myCluster = parcluster('local');
% for Name = {'hiv','threestep'}
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
