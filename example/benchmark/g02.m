function [f, g] = g02(x)

g(1) = 1;
g(2) = 0;
f1 = 0;
f2 = 1;
f3 = 0;
for i=1:length(x)
    g(1) = g(1) * x(i);
    g(2) = g(2) + x(i);
    f1 = f1 + cos(x(i)) ^ 4;
    f2 = f2 * cos(x(i)) ^ 2;
    f3 = f3 + i * x(i) * x(i);
end

g(1) = 0.75 - g(1);
g(2) = g(2) - 7.5 * length(x);
f = - abs( ( f1 - 2 * f2 ) / sqrt(f3) );