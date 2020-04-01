function f = k_tablet(x)
% Unconstrained benchmark function k_tablet.
% 
% The optimum is located at x* = (0, ..., 0) where f(x*) = 0.
% 
% [SYNTAX]
% [f, g] = k_tablet(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


n = length(x);
k = 0.25 * n;

f = sum( x(1:k) .^ 2 ) + sum( ( 100 * x(k+1:n) ) .^ 2 );
