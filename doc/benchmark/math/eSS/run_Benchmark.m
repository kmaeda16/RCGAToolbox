
n_repeat = 10;

for i = 1 : n_repeat
%     Benchmark_eSS(i); % For normal calculation
    batch(@Benchmark_eSS,0,{i}); % For batch calculation
end
