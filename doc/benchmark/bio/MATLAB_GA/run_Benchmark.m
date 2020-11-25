
n_repeat = 5;
dirname = 'Results';

mkdir(dirname);
cd(dirname);
addpath('..');


%% Start Calculation
for Name = {'hiv','threestep'}
    name = char(Name);
    for i = 1 : n_repeat
%         Biological_MATLAB_GA(name,i); % For normal calculation
        batch(@Biological_MATLAB_GA,0,{name,i}); % For batch calculation
    end
end

rmpath('..');
cd('..');


%% Save diary
% jobid = 30;
% for Name = {'hiv','threestep'}
%     name = char(Name);
%     for j = 1 : n_repeat
%         jobid = jobid + 1;
%         jobid = sprintf('job%d',jobid);
%         diaryname = sprintf('MATLAB_GA_%s_diary_%d.txt',name,j);
%         expression = sprintf('diary(%s,''%s'')',jobid,diaryname);
%         eval(expression);
%     end
% end
