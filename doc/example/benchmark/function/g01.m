function [f, g] = g01(x)

f = 5 * sum( x(1:4) ) - 5 * sum( x(1:4) .^2 ) - sum( x(5:13) );

g(1) = 2 * x(1) + 2 * x(2) +  x(10) + x(11) - 10;
g(2) = 2 * x(1) + 2 * x(3) +  x(10) + x(12) - 10;
g(3) = 2 * x(2) + 2 * x(3) + x(11) + x(12) - 10;
g(4) = - 8 * x(1) + x(10);
g(5) = - 8 * x(2) + x(11);
g(6) = - 8 * x(3) + x(12);
g(7) = - 2 * x(4) - x(5) + x(10);
g(8) = - 2 * x(6) - x(7) + x(11);
g(9) = - 2 * x(8) - x(9) + x(12);