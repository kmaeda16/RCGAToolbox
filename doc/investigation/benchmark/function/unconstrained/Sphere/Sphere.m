function f = Sphere(x)
% Unconstrained benchmark function Sphere.
% 
% The optimum is located at x* = (0, ..., 0) where f(x*) = 0.
% 
% [SYNTAX]
% [f, g] = Sphere(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


f = sum(x.^2);
