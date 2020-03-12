function f = MMbenchmark(x)
% Unconstrained benchmark function MMbenchmark.
% 
% The optimum is located at x* = (0.699, ..., 0.699) where f(x*) = 0.
% 
% [SYNTAX]
% [f, g] = MMbenchmark(x)
% 
% [INPUT]
% x : Decision variables
% 
% [OUTPUT]
% f : Objective function


n = length(x);
f = 0;
for i = 1 : 0.5 * n
    f = f ...
        + ( 10^x(2*i-1) * 0.1 / ( 10^(x(2*i)-2.0) + 0.1 ) - 3.33333333333333 ) ^ 2 ...
        + ( 10^x(2*i-1) * 5.0 / ( 10^(x(2*i)-2.0) + 5.0 ) - 4.95049504950495 ) ^ 2;
end
