
n_repeat = 10;

for i = 1 : n_repeat
%     doBenchmark_eSS(i); % For normal calculation
    batch(@doBenchmark_eSS,0,{i}); % For batch calculation
end
