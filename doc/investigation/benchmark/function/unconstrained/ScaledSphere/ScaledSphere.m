function f = ScaledSphere(x)
% Unconstrained benchmark function ScaledSphere.
% 
% The optimum is located at x* = (0, ..., 0) where f(x*) = 0.
% 
% [SYNTAX]
% [f, g] = ScaledSphere(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


n = length(x);
i = 1 : n;
f = sum( ( i .* x ) .^ 2 );
