function [f, g] = g11(x)

DELTA = 1e-4;

f = x(1) ^ 2 + ( x(2) - 1 ) ^ 2;

g(1) = abs( x(2) - x(1) ^ 2 ) - DELTA;
