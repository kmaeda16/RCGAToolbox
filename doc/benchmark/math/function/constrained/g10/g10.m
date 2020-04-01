function [f, g] = g10(x)
% Constrained benchmark function g10.
% 
% The optimum is located at x* = (579.3167, 1359.943, 5110.071, 182.0174,
% 295.5985, 217.9799, 286.4162, 395.5979) where f(x*) = 7049.3307.
% 
% [SYNTAX]
% [f, g] = g10(x)
% 
% [INPUT]
% x : Decision variables (8 dimensional)
% 
% [OUTPUT]
% f : Objective function
% g : Constraint functions (6 dimensional)


f = x(1) + x(2) + x(3);

g(1) = - 1 + 0.0025 * ( x(4) + x(6) );
g(2) = - 1 + 0.0025 * ( x(5) + x(7) - x(4) );
g(3) = - 1 + 0.01 * ( x(8) - x(5) );
g(4) = - x(1) * x(6) + 833.33252 * x(4) + 100 * x(1) - 83333.333;
g(5) = - x(2) * x(7) + 1250 * x(5) + x(2) * x(4) - 1250 * x(4);
g(6) = - x(3) * x(8) + 1250000 + x(3) * x(5) - 2500 * x(5);
