function [f, g] = g13(x)

DELTA = 1e-4;

f = exp( prod(x) );

g(1) = abs( sum( x .^ 2 ) - 10 ) - DELTA;
g(2) = abs( x(2) * x(3) - 5 * x(4) * x(5) ) - DELTA;
g(3) = abs( x(1) ^ 3 + x(2) ^ 3 + 1 ) - DELTA;
