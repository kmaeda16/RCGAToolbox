
n_repeat = 10;

for name = {'hiv','threestep'}
    for i = 1 : n_repeat
%         doBiological_RCGA(name(i),i); % For normal calculation
        batch(@doBiological_RCGA,0,{name(i),i}); % For batch calculation
    end
end
