function f = Rastrigin(x)

n = length(x);
f = 10 * n;
for i = 1 : n
    f = f + x(i)^2 - 10 * cos( 2 * pi * x(i) );
end
