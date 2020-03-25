function [f, g] = g08(x)
% Constrained benchmark function g08.
% 
% The optimum is located at x* = (1.2279713, 4.2453733) where f(x*) =
% -0.095825.
% 
% [SYNTAX]
% [f, g] = g08(x)
% 
% [INPUT]
% x : Decision variables (2 dimensional)
% 
% [OUTPUT]
% f : Objective function
% g : Constraint functions (2 dimensional)


f = - sin( 2 * pi * x(1) ) ^ 3 * sin( 2 * pi * x(2) ) / ( x(1) ^ 3 * ( x(1) + x(2) ) );

g(1) = x(1) ^ 2 - x(2) + 1;
g(2) = 1 - x(1) + ( x(2) - 4 ) ^ 2 ;
