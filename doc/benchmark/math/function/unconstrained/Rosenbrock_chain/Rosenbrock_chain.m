function f = Rosenbrock_chain(x)
% Unconstrained benchmark function Rosenbrock_chain.
% 
% The optimum is located at x* = (1, ..., 1) where f(x*) = 0.
% 
% [SYNTAX]
% [f, g] = Rosenbrock_chain(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


n = length(x);
f = 0;
for i = 1 : n - 1
    f = f + 100 * ( x(i+1) - x(i)^2 ) ^ 2 + ( 1 - x(i) ) ^ 2;
end
