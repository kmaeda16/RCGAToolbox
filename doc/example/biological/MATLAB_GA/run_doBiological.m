
n_repeat = 5;

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

% diaryname = sprintf('MATLAB_GA_%s_diary_%d.dat',problem_name,idum);
% diary(job433,diaryname);
