function [f, g] = fitnessExample(x)
% fitnessExample is an example of fitness function.
% 
% [SYNTAX]
% [f, g] = fitnessExample(x)
% 
% [INPUT]
% x : Decision variables.
% 
% [OUTPUT]
% f : Objective function to be miniimzed.
% g : Constraint functions to be less than zero.
% 
% --------------------- Example Problem ---------------------
% Minimize:
%   f = x(1)^2 + x(2)^2 + ... + x(10)^2
% 	
% Subject to:
%   g(1) = x(1) * x(2) + 1 <= 0
%   g(2) = x(1) + x(2) + 1 <= 0
% 	-5.12 <= x(i) <= 5.12 for all i
% -----------------------------------------------------------


f    = sum( x .^ 2 );

g(1) = x(1) * x(2) + 1;
g(2) = x(1) + x(2) + 1;
