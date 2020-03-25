
n_repeat = 5;

for Name = {'hiv','threestep'}
    name = char(Name);
    for i = 1 : n_repeat
        dirname = sprintf('result_%s_%d',name,i);
        mkdir(dirname);
        cd(dirname);
        addpath('..');
%         doBiological_eSS(name,i); % For normal calculation
        batch(@doBiological_eSS,0,{name,i}); % For batch calculation
        rmpath('..');
        cd('..');
    end
end
