function f = Ackley(x)

n = length(x);

temp1 = 0;
temp2 = 0;
for i = 1 : n
    temp1 = temp1 + x(i) * x(i);
    temp2 = temp2 + cos( 2.0 * pi * x(i) );
end
f = 20.0 - 20.0 * exp( - 0.2 * sqrt( 1.0 / n * temp1 ) ) ...
    + exp(1) - exp( 1.0 / n * temp2);
