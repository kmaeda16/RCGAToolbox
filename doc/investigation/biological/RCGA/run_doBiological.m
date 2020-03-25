
n_repeat = 5;

for Name = {'hiv','threestep'}
    name = char(Name);
    for i = 1 : n_repeat
        dirname = 'Results';
        mkdir(dirname);
        cd(dirname);
        addpath('..');
%         doBiological_RCGA_UNDXMGG(name,i); % For normal calculation
        batch(@doBiological_RCGA_UNDXMGG,0,{name,i}); % For batch calculation
%         doBiological_RCGA_REXstarJGG(name,i); % For normal calculation
        batch(@doBiological_RCGA_REXstarJGG,0,{name,i}); % For batch calculation
        rmpath('..');
        cd('..');
    end
end
