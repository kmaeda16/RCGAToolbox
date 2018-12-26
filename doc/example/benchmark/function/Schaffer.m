function f = Schaffer(x)

n = length(x);
f = 0;
for i = 1 : n - 1
    temp1 = ( x(i)^2 + x(i+1)^2 ) ^ 0.25 ;
    temp2 = 50.0 * ( x(i)^2 + x(i+1)^2 ) ^ 0.1 ;
    f = f + temp1 * ( sin(temp2)^2 + 1.0 );
end
