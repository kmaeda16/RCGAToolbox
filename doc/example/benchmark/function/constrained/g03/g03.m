function [f, g] = g03(x)
% Constrained benchmark function g03.
% 
% The optimum is located at x_i* = 1 / sqrt(n) (i = 1,... , n) where f(x*)
% = -1
% 
% [SYNTAX]
% [f, g] = g03(x)
% 
% [INPUT]
% x : Decision variables (10 dimensional)
% 
% [OUTPUT]
% f : Objective function
% g : Constraint functions (1 dimensional)


n = length(x);

DELTA = 1e-4;

f = - sqrt(n) ^ n * prod(x);

g(1) = abs( sum( x .^ 2 ) - 1 ) - DELTA;
