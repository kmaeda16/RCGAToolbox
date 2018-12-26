function f = Rosenbrock_star(x)

n = length(x);
f = 0;
for i = 2 : n
    f = f + 100.0 * ( x(1) - x(i)^2 ) ^ 2 + ( 1.0 - x(i) ) ^ 2;
end
