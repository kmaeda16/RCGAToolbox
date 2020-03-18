function [c, ceq] = cst_wrapper_nodecoding(fitnessfun, x)
% Constrained benchmark function g01.
% 
% The optimum is located at x* = (1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 1)
% where f(x*) = -15.
% 
% [SYNTAX]
% [f, g] = g01(x)
% 
% [INPUT]
% x : Decision variables (13 dimensional)
% 
% [OUTPUT]
% f : Objective function
% g : Constraint functions (9 dimensional)

% x = decodingfun(gene);
[~, c] = fitnessfun(x);
ceq = [];
