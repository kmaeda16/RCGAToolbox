function [f, g] = g13(x)
% Constrained benchmark function g13.
% 
% The optimum is located at x* = (-1.717143, 1.595709, 1.827247,
% -0.7636413, -0.763645) where f(x*) = 0.0539498.
% 
% [SYNTAX]
% [f, g] = g13(x)
% 
% [INPUT]
% x : Decision variables (5 dimensional)
% 
% [OUTPUT]
% f : Objective function
% g : Constraint functions (3 dimensional)


DELTA = 1e-4;

f = exp( prod(x) );

g(1) = abs( sum( x .^ 2 ) - 10 ) - DELTA;
g(2) = abs( x(2) * x(3) - 5 * x(4) * x(5) ) - DELTA;
g(3) = abs( x(1) ^ 3 + x(2) ^ 3 + 1 ) - DELTA;
