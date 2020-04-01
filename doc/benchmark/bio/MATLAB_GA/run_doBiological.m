
n_repeat = 5;

%%
for Name = {'hiv','threestep'}
    name = char(Name);
    for i = 1 : n_repeat
        dirname = 'Results';
        mkdir(dirname);
        cd(dirname);
        addpath('..');
%         doBiological_MATLAB_GA(name,i); % For normal calculation
        batch(@doBiological_MATLAB_GA,0,{name,i}); % For batch calculation
        rmpath('..');
        cd('..');
    end
end

%%
i = 30;
for Name = {'hiv','threestep'}
    name = char(Name);
    for j = 1 : n_repeat
        i = i + 1;
        jobid = sprintf('job%d',i);
        diaryname = sprintf('MATLAB_GA_%s_diary_%d.txt',name,j);
        expression = sprintf('diary(%s,''%s'')',jobid,diaryname);
        eval(expression);
    end
end
