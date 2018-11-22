function f = ScaledSphere(x)

n = length(x);
f = 0;
for i = 1 : n
    f = f + ( i * x(i) ) ^ 2 ;
end
