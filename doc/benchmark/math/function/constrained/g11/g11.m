function [f, g] = g11(x)
% Constrained benchmark function g11.
% 
% The optimum is located at x* = (+-1/sqrt(2), 1/2) where f(x*) = 0.75.
% 
% [SYNTAX]
% [f, g] = g11(x)
% 
% [INPUT]
% x : Decision variables (2 dimensional)
% 
% [OUTPUT]
% f : Objective function
% g : Constraint functions (1 dimensional)


DELTA = 1e-4;

f = x(1) ^ 2 + ( x(2) - 1 ) ^ 2;

g(1) = abs( x(2) - x(1) ^ 2 ) - DELTA;
