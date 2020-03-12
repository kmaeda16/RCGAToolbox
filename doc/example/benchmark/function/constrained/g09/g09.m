function [f, g] = g09(x)
% Constrained benchmark function g09.
% 
% The optimum is located at x* = (2.330499, 1.951372, -0.4775414, 4.365726,
% -0.6244870, 1.038131, 1.594227) where f(x*) = 680.6300573.
% 
% [SYNTAX]
% [f, g] = g09(x)
% 
% [INPUT]
% x : Decision variables (7 dimensional)
% 
% [OUTPUT]
% f : Objective function
% g : Constraint functions (4 dimensional)


f = ( x(1) - 10 ) ^ 2 + 5 * ( x(2) - 12 ) ^ 2 + x(3) ^ 4 + 3 * ( x(4) - 11 ) ^ 2 ...
    + 10 * x(5) ^ 6 + 7 * x(6) ^ 2 + x(7) ^ 4 - 4 * x(6) * x(7) - 10 * x(6) - 8 * x(7);

g(1) = - 127 + 2 * x(1) ^ 2 + 3 * x(2) ^ 4 + x(3) + 4 * x(4) ^ 2 + 5 * x(5);
g(2) = - 282 + 7 * x(1) + 3 * x(2) + 10 * x(3) ^ 2 + x(4) - x(5);
g(3) = - 196 + 23 * x(1) + x(2) ^ 2 + 6 * x(6) ^ 2 - 8 * x(7);
g(4) = 4 * x(1) ^ 2 + x(2) ^ 2 - 3 * x(1) * x(2) + 2 * x(3) ^ 2 + 5 * x(6) - 11 * x(7);
