function f = Rosenbrock_star(x)
% Unconstrained benchmark function Rosenbrock_star.
% 
% The optimum is located at x* = (1, ..., 1) where f(x*) = 0.
% 
% [SYNTAX]
% [f, g] = Rosenbrock_star(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


n = length(x);
f = 0;
for i = 2 : n
    f = f + 100.0 * ( x(1) - x(i)^2 ) ^ 2 + ( 1.0 - x(i) ) ^ 2;
end
