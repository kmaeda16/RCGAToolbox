
n_repeat = 1;

for i = 1 : n_repeat
    doBenchmark_RCGA(i); % For normal calculation
%     batch(@doBenchmark,0,{i}); % For batch calculation
end
