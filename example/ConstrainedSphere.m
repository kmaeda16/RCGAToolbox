function [f, g] = ConstrainedSphere(x)

f = sum( x .^ 2 );
g(1) = x(1) * x(2) + 1;
g(2) = x(1) + x(2) + 1;
