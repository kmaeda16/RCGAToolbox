function f = Rastrigin(x)
% Unconstrained benchmark function Rastrigin.
% 
% The optimum is located at x* = (0, ..., 0) where f(x*) = 0.
% 
% [SYNTAX]
% [f, g] = Rastrigin(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


n = length(x);
f = 10 * n;
for i = 1 : n
    f = f + x(i) ^ 2 - 10 * cos( 2 * pi * x(i) );
end
