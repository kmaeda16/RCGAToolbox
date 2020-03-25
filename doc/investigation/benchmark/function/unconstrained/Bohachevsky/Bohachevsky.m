function f = Bohachevsky(x)
% Unconstrained benchmark function Bohachevsky.
% 
% The optimum is located at x* = (0, ..., 0) where f(x*) = 0.
% 
% [SYNTAX]
% [f, g] = Bohachevsky(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


n = length(x);

f = 0;
for i = 1 : n - 1;
    f = f + x(i)^2 + 2.0 * x(i+1)^2 - 0.3 * cos( 3.0 * pi * x(i) ) ...
        - 0.4 * cos( 4.0 * pi * x(i+1) ) + 0.7;
end
