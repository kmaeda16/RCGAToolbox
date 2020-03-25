function f = Schwefel(x)
% Unconstrained benchmark function Schwefel.
% 
% The optimum is located at x* = (-420.9687, ..., -420.9687) where f(x*) =
% 0.
% 
% [SYNTAX]
% [f, g] = Schwefel(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


n = length(x);
f = 418.9828873 * n;
for i = 1 : n
    f = f + x(i) * sin( sqrt( abs( x(i) ) ) );
end
