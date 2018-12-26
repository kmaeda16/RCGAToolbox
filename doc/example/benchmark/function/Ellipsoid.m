function f = Ellipsoid(x)

n = length(x);
f = 0;
for i = 1 : n
    f = f + ( 1000 ^ ( ( i - 1 ) / ( n - 1 ) ) * x(i) ) ^ 2;
end