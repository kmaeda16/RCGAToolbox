function f = Cigar(x)
% Unconstrained benchmark function Cigar.
% 
% The optimum is located at x* = (0, ..., 0) where f(x*) = 0.
% 
% [SYNTAX]
% [f, g] = Cigar(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


f = x(1) ^ 2 + sum( ( 1000 * x(2:end) ) .^ 2 );
