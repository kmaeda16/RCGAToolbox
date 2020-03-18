
n_repeat = 10;

for i = 1 : n_repeat
%     doBenchmark_RCGA(i); % For normal calculation
    batch(@doBenchmark_RCGA,0,{i}); % For batch calculation
end
