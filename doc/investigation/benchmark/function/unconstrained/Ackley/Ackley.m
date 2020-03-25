function f = Ackley(x)
% Unconstrained benchmark function Ackley.
% 
% The optimum is located at x* = (0, ..., 0) where f(x*) = 0.
% 
% [SYNTAX]
% [f, g] = Ackley(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


n = length(x);

temp1 = 0;
temp2 = 0;
for i = 1 : n
    temp1 = temp1 + x(i) * x(i);
    temp2 = temp2 + cos( 2.0 * pi * x(i) );
end
f = 20.0 - 20.0 * exp( - 0.2 * sqrt( 1.0 / n * temp1 ) ) ...
    + exp(1) - exp( 1.0 / n * temp2 );
