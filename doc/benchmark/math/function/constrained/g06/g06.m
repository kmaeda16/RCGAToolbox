function [f, g] = g06(x)
% Constrained benchmark function g06.
% 
% The optimum is located at x* = (14.095, 0.84296) where f(x*) =
% -6961.81388.
% 
% [SYNTAX]
% [f, g] = g06(x)
% 
% [INPUT]
% x : Decision variables (2 dimensional)
% 
% [OUTPUT]
% f : Objective function
% g : Constraint functions (2 dimensional)


f = ( x(1) - 10 ) ^ 3 + ( x(2) - 20 ) ^ 3;

g(1) = - ( x(1) - 5 ) ^ 2 - ( x(2) - 5 ) ^ 2 + 100;
g(2) =   ( x(1) - 6 ) ^ 2 + ( x(2) - 5 ) ^ 2 - 82.81;
