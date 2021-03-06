function [f, g] = g05(x)
% Constrained benchmark function g05.
% 
% The known best value is x* = (679.9453, 1026.067, 0.1188764, -0.3962336)
% where f(x*) = 5126.4981.
% 
% [SYNTAX]
% [f, g] = g05(x)
% 
% [INPUT]
% x : Decision variables (4 dimensional)
% 
% [OUTPUT]
% f : Objective function
% g : Constraint functions (5 dimensional)


DELTA = 1e-4;

f = 3 * x(1) + 0.000001 * x(1) ^ 3 + 2 * x(2) + ( 0.000002 / 3 ) * x(2) ^ 3;

g(1) = - x(4) + x(3) - 0.55;
g(2) = - x(3) + x(4) - 0.55;
g(3) = abs( 1000 * sin( - x(3) - 0.25 ) + 1000 * sin( - x(4)        - 0.25 ) + 894.8 - x(1) ) - DELTA;
g(4) = abs( 1000 * sin(   x(3) - 0.25 ) + 1000 * sin(   x(3) - x(4) - 0.25 ) + 894.8 - x(2) ) - DELTA;
g(5) = abs( 1000 * sin(   x(4) - 0.25 ) + 1000 * sin(   x(4) - x(3) - 0.25 ) + 1294.8 ) - DELTA;
