function [f, g] = g02(x)
% Constrained benchmark function g02.
% 
% The optimum is unknown. The known best value is f(x*) = -0.803619.
% 
% [SYNTAX]
% [f, g] = g02(x)
% 
% [INPUT]
% x : Decision variables (20 dimensional)
% 
% [OUTPUT]
% f : Objective function
% g : Constraint functions (2 dimensional)


g(1) = 1;
g(2) = 0;
f1 = 0;
f2 = 1;
f3 = 0;

for i = 1 : length(x)
    g(1) = g(1) * x(i);
    g(2) = g(2) + x(i);
    f1 = f1 + cos( x(i) ) ^ 4;
    f2 = f2 * cos( x(i) ) ^ 2;
    f3 = f3 + i * x(i) * x(i);
end

g(1) = 0.75 - g(1);
g(2) = g(2) - 7.5 * length(x);
f = - abs( ( f1 - 2 * f2 ) / sqrt(f3) );
