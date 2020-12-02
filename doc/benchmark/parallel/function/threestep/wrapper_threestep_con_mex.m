function [f, g] = wrapper_threestep_con_mex(x)
% Solution:
% x = log10([1, 1, 2, 1, 2, 1, 1, 1, 2, 1, 2, 1, 1, 1, 2, 1, 2, 1, 0.1, ...
% 1, 0.1, 0.1, 1, 0.1, 0.1, 1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]);

x = 10 .^ x;

[f, g] = threestep_con_mex(x);
