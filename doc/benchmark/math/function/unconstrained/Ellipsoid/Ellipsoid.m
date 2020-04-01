function f = Ellipsoid(x)
% Unconstrained benchmark function Ellipsoid.
% 
% The optimum is located at x* = (0, ..., 0) where f(x*) = 0.
% 
% [SYNTAX]
% [f, g] = Ellipsoid(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


n = length(x);
i = 1 : n;
f = sum( ( 1000 .^ ( ( i - 1 ) ./ ( n - 1 ) ) .* x ) .^ 2 );
