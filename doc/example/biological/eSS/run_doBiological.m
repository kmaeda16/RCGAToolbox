
n_repeat = 10;

for name = {'hiv','threestep'}
    for i = 1 : n_repeat
%         doBiological_eSS(name(i),i); % For normal calculation
        batch(@doBiological_eSS,0,{name(i),i}); % For batch calculation
    end
end
