
n_repeat = 10;

for i = 1 : n_repeat
%     doBenchmark_MATLAB_GA(i); % For normal calculation
    batch(@doBenchmark_MATLAB_GA,0,{i}); % For batch calculation
end
