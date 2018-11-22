function [f, g] = g03(x)

n = length(x);

DELTA = 1e-4;

f = - sqrt(n) ^ n * prod(x);

g(1) = abs( sum( x .^ 2 ) - 1 ) - DELTA;
